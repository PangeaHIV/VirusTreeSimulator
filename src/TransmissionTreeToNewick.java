import dr.app.tools.NexusExporter;
import dr.app.util.Arguments;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.ExponentialGrowth;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.epidemiology.LogisticGrowthN0;
import dr.math.MathUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Output a CSV transmission tree as Newick format
 *
 * @author mhall
 */

public class TransmissionTreeToNewick {


    public static final String IDREC = "Host_ID";
    public static final String IDTR = "Infector_ID";
    public static final String TIME_INF = "Infection_Time";
    public static final String TIME_REC = "Removal_Time";

    private ArrayList<InfectedUnit> units;
    private HashMap<String, InfectedUnit> idMap;
    private String outputFileRoot;


    public TransmissionTreeToNewick(String fileName, String outputFileRoot){
        units = new ArrayList<>();
        idMap = new HashMap<>();
        this.outputFileRoot = outputFileRoot;
        try {
            readEvents(fileName);
        } catch(IOException e){
            e.printStackTrace();
        }
    }

    private enum EventType{
        INFECTION, SAMPLE
    }


    private void run() throws IOException{
        FlexibleTree tree = makeTree();

        NexusExporter exporter = new NexusExporter(new PrintStream(outputFileRoot + ".nex"));
        exporter.exportTree(tree);
    }

    private void readEvents(String fileName) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(fileName));

        String[] headers = reader.readLine().split(",");

        int infecteeColumn = -1;
        int infectorColumn = -1;
        int infectionTimeColumn = -1;
        int endTimeColumn = -1;

        for(int i=0; i<headers.length; i++){
            String headerItem = headers[i].replaceAll("\"", "");
            if(headerItem.equals(IDREC)){
                infecteeColumn = i;
            } else if(headerItem.equals(IDTR)){
                infectorColumn = i;
            } else if(headerItem.equals(TIME_INF)){
                infectionTimeColumn = i;
            } else if(headerItem.equals(TIME_REC)){
                endTimeColumn = i;
            }
        }

        if(infecteeColumn == -1 || infectorColumn == -1 || infectionTimeColumn == -1 || endTimeColumn == -1){
            throw new RuntimeException("Not all required columns are present in the file");
        }

        String line = reader.readLine();

        while(line!=null){
            String[] entries = line.split(",");

            InfectedUnit infectee = new InfectedUnit("ID_"+entries[infecteeColumn]);

            units.add(infectee);
            idMap.put("ID_"+entries[infecteeColumn], infectee);

            line = reader.readLine();
        }

        reader = new BufferedReader(new FileReader(fileName));

        reader.readLine();
        line = reader.readLine();

        while(line!=null){
            String[] entries = line.split(",");

            InfectedUnit infectee = idMap.get("ID_"+entries[infecteeColumn]);

            if(!idMap.containsKey("ID_"+entries[infectorColumn]) & !entries[infectorColumn].equals("NA")){
                throw new RuntimeException(entries[infectorColumn] + " does not appear in the infectee column of "
                        + fileName);
            }

            if(!entries[infectorColumn].equals("NA")) {
                InfectedUnit infector = idMap.get("ID_" + entries[infectorColumn]);

                Event infection = new Event(EventType.INFECTION, Double.parseDouble(entries[infectionTimeColumn]), infector,
                        infectee);

                infector.addInfectionEvent(infection);
                infectee.setInfectionEvent(infection);

                infectee.parent = infector;
            } else {
                Event infection = new Event(EventType.INFECTION, Double.parseDouble(entries[infectionTimeColumn]), null, infectee);
                infectee.setInfectionEvent(infection);
            }

            infectee.addSamplingEvent(Double.parseDouble(entries[endTimeColumn]), 1);

            line = reader.readLine();
        }

    }

    private FlexibleTree makeTree(){

        // find the first case

        ArrayList<InfectedUnit> introducedCases = new ArrayList<>();

        for(InfectedUnit unit : units){
            if(unit.parent==null){
                introducedCases.add(unit);
            }
        }

        if(introducedCases.size()==0){
            throw new RuntimeException("Can't find a first case");
        }

        if(introducedCases.size()>1){
            throw new RuntimeException("We require a single, connected tree");
        }

        double earliestEventTime = Double.POSITIVE_INFINITY;
        double latestEventTime = Double.NEGATIVE_INFINITY;

        for(InfectedUnit unit : units){
            if(unit.infectionEvent.time < earliestEventTime){
                earliestEventTime = unit.infectionEvent.time;
            }
            for(Event child : unit.childEvents){
                if(child.time > latestEventTime){
                    latestEventTime = child.time;
                }
            }
        }

        HashMap<Event, SimpleNode> eventToNode = new HashMap<>();

        ArrayList<Event> allInfections = new ArrayList<>();

        for(InfectedUnit unit : units){
            if(!(allInfections.contains(unit.infectionEvent))) {
                allInfections.add(unit.infectionEvent);
            }
        }

        Collections.sort(allInfections);

        SimpleTree out = null;

        for(Event infectionEvent: allInfections){
            InfectedUnit infectee = infectionEvent.infectee;
            ArrayList<Event> infecteeChildEvents = infectee.childEvents;
            Collections.sort(infecteeChildEvents);
            SimpleNode eventNode;
            if(infectionEvent == allInfections.get(0)){
                double nodeTime = latestEventTime - infectionEvent.time;
                eventNode = new SimpleNode();
                out = new SimpleTree(eventNode);
                out.beginTreeEdit();
                eventNode.setHeight(nodeTime);
                eventToNode.put(infectionEvent, eventNode);
            } else {
                if(!(eventToNode.containsKey(infectionEvent))){
                    throw new RuntimeException("No record of infection node that should already have been generated.");
                }
                eventNode =  eventToNode.get(infectionEvent);
            }
            SimpleNode lastEventNode = eventNode;

            for(Event childEvent : infecteeChildEvents){
                double nodeTime = latestEventTime - childEvent.time;
                eventNode = new SimpleNode();
                eventToNode.put(childEvent, eventNode);
                out.addChild(lastEventNode, eventNode);
                eventNode.setHeight(nodeTime);
                lastEventNode = eventNode;
                if(childEvent.type == EventType.SAMPLE){
                    eventNode.setTaxon(new Taxon(infectee.id));
                }

            }

        }
        out.endTreeEdit();

        return new FlexibleTree(new FlexibleNode(out, out.getRoot(), true));
    }






    private class InfectedUnit{
        private String id;
        private ArrayList<Event> childEvents;
        private Event infectionEvent;
        private InfectedUnit parent;

        private InfectedUnit(String id){
            this.id = id;
            parent = null;
            childEvents = new ArrayList<>();
        }

        private void addSamplingEvent(double time){
            if(infectionEvent!=null && time < infectionEvent.time){
                throw new RuntimeException("Adding an event to case "+id+" before its infection time");
            }
            childEvents.add(new Event(EventType.SAMPLE, time));
        }

        private void addSamplingEvent(double time, int instances){
            if(infectionEvent!=null && time < infectionEvent.time){
                throw new RuntimeException("Adding an event to case "+id+" before its infection time");
            }
            childEvents.add(new Event(EventType.SAMPLE, time, instances));
        }

        private void setInfectionEvent(double time, InfectedUnit infector){
            setInfectionEvent(new Event(EventType.INFECTION, time, infector, this));
        }

        private void setInfectionEvent(Event event){
            for(Event childEvent : childEvents){
                if(event.time > childEvent.time){

                    if(childEvent.type == EventType.SAMPLE){
                        throw new RuntimeException("Setting infection time for case "+id+" after its sampling at "+
                                childEvent.time);
                    } else {
                        String childUnitName = childEvent.infectee.id;

                        throw new RuntimeException("Setting infection time for case "+id+" after it infected "
                                +childUnitName+" at "+childEvent.time);
                    }
                }
            }
            infectionEvent = event;
        }

        private void addChildInfectionEvent(double time, InfectedUnit infectee){
            addInfectionEvent(new Event(EventType.INFECTION, time, this, infectee));
        }

        private void addInfectionEvent(Event event){
            if(infectionEvent!=null && event.time < infectionEvent.time){
                throw new RuntimeException("Adding an infection event to case "+id+" at "+event.time+" before its " +
                        "infection time at "+infectionEvent.time);
            }
            childEvents.add(event);
        }

        private void sortEvents(){
            Collections.sort(childEvents);
            Collections.reverse(childEvents);
        }

    }

    private class Event implements Comparable<Event>{

        private EventType type;
        private double time;
        // multiple samples at the same time and multiple transmissions to the same host at the same time are now _one_ event
        private int instances;
        private InfectedUnit infector;
        private InfectedUnit infectee;


        private Event(EventType type, double time){
            this.type = type;
            this.time = time;
            this.instances = 1;
        }

        private Event(EventType type, double time, InfectedUnit infector, InfectedUnit infectee){
            this.type = type;
            this.time = time;
            this.infector = infector;
            this.infectee = infectee;
            this.instances = 1;
        }

        private Event(EventType type, double time, int instances){
            this.type = type;
            this.time = time;
            this.instances = instances;
        }

        private Event(EventType type, double time, int instances, InfectedUnit infector, InfectedUnit infectee){
            this.type = type;
            this.time = time;
            this.infector = infector;
            this.infectee = infectee;
            this.instances = instances;
        }

        public void addInstance(){
            this.instances++;
        }

        public void setInstanceCount(int count){
            this.instances = count;
        }


        public int compareTo(Event event) {
            return Double.compare(time, event.time);
        }
    }


    public static void main(String[] args){


        String infectionsFileName = args[0];
        String outputFileRoot = args[1];


        TransmissionTreeToNewick instance = new TransmissionTreeToNewick(infectionsFileName, outputFileRoot);

        try {
            instance.run();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
