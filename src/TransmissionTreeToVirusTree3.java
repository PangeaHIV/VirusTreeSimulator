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
 * Simulated a virus tree given a transmission tree and dates of sampling
 *
 * @author mhall
 */

public class TransmissionTreeToVirusTree3 {

    protected static PrintStream progressStream = System.out;

    private enum ModelType{CONSTANT, EXPONENTIAL, LOGISTIC}

    public static final String HELP = "help";

    public static final String DEMOGRAPHIC_MODEL = "demoModel";

    public static final String[] demographics = {"Constant", "Exponential", "Logistic"};

    public static final String STARTING_POPULATION_SIZE = "N0";
    public static final String GROWTH_RATE = "growthRate";
    public static final String T50 = "t50";
    public static final String SEED = "seed";
    public static final String FORCE_COALESCENCE = "forceCoalescence";


    public static final String IDREC = "IDREC";
    public static final String IDTR = "IDTR";
    public static final String TIME_TR = "TIME_TR";
    public static final String IDPOP = "IDPOP";
    public static final String TIME_SEQ = "TIME_SEQ";
    public static final String SEQ_COUNT = "SEQ_COUNT";


    private DemographicFunction demFunct;
    private ArrayList<InfectedUnit> units;
    private HashMap<String, InfectedUnit> idMap;
    private String outputFileRoot;

    private double coalescentProbability;

    public TransmissionTreeToVirusTree3(String fileName,
                                        DemographicFunction demFunct, String outputFileRoot){
        this.demFunct = demFunct;
        units = new ArrayList<InfectedUnit>();
        idMap = new HashMap<String, InfectedUnit>();
        this.outputFileRoot = outputFileRoot;
        coalescentProbability = 1;
        try {
            readSamplingEvents(fileName);
            readInfectionEvents(fileName);

        } catch(IOException e){
            e.printStackTrace();
        }
    }

    public TransmissionTreeToVirusTree3(String sampFileName, String transFileName,
                                        DemographicFunction demFunct, String outputFileRoot){
        this.demFunct = demFunct;
        units = new ArrayList<InfectedUnit>();
        idMap = new HashMap<String, InfectedUnit>();
        this.outputFileRoot = outputFileRoot;
        try {
            readInfectionEvents(transFileName);
            readSamplingEvents(sampFileName);

        } catch(IOException e){
            e.printStackTrace();
        }
    }



    private enum EventType{
        INFECTION, SAMPLE
    }


    private void run(boolean forceCoalescences) throws IOException{
        ArrayList<FlexibleTree> detailedTrees = makeTrees(forceCoalescences);
        ArrayList<FlexibleTree> simpleTrees = new ArrayList<FlexibleTree>();

        for(FlexibleTree tree : detailedTrees) {
            FlexibleTree wbTree = makeWellBehavedTree(tree);
            wbTree.setAttribute("firstCase", tree.getAttribute("firstCase"));

            simpleTrees.add(wbTree);
        }


        for(FlexibleTree tree: detailedTrees){
            NexusExporter exporter = new NexusExporter(new PrintStream(outputFileRoot
                    + tree.getAttribute("firstCase") + "_detailed.nex"));
            exporter.exportTree(tree);
        }

        for(FlexibleTree tree: simpleTrees){
            NexusExporter exporter = new NexusExporter(new PrintStream(outputFileRoot
                    + tree.getAttribute("firstCase") + "_simple.nex"));
            exporter.exportTree(tree);
        }

    }

    private void readInfectionEvents(String fileName) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(fileName));

        String[] headers = reader.readLine().split(",");

        int infecteeColumn = -1;
        int infectorColumn = -1;
        int timeColumn = -1;

        for(int i=0; i<headers.length; i++){
            String headerItem = headers[i].replaceAll("\"", "");
            if(headerItem.equals(IDREC)){
                infecteeColumn = i;
            } else if(headerItem.equals(IDTR)){
                infectorColumn = i;
            } else if(headerItem.equals(TIME_TR)){
                timeColumn = i;
            }
        }

        if(infecteeColumn == -1 || infectorColumn == -1 || timeColumn == -1){
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
                throw new RuntimeException(entries[infectorColumn] + "does not appear in the infectee column of "
                        + fileName);
            }

            if(!entries[infectorColumn].equals("NA")) {
                InfectedUnit infector = idMap.get("ID_" + entries[infectorColumn]);

                Event infection = new Event(EventType.INFECTION, Double.parseDouble(entries[timeColumn]), infector,
                        infectee);

                infector.addInfectionEvent(infection);
                infectee.setInfectionEvent(infection);

                infectee.parent = infector;
            } else {
                Event infection = new Event(EventType.INFECTION, Double.parseDouble(entries[timeColumn]), null, infectee);
                infectee.setInfectionEvent(infection);
            }

            line = reader.readLine();
        }

    }

    private void readSamplingEvents(String fileName) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(fileName));

        String[] headers = reader.readLine().split(",");

        int unitColumn = -1;
        int samplingTimeColumn = -1;
        int sampleCountColumn = -1;

        for(int i=0; i<headers.length; i++){
            String headerItem = headers[i].replaceAll("\"", "");
            if(headerItem.equals(IDPOP)){
                unitColumn = i;
            } else if(headerItem.equals(TIME_SEQ)){
                samplingTimeColumn = i;
            } else if(headerItem.equals(SEQ_COUNT)){
                sampleCountColumn = i;
            }
        }

        if(unitColumn == -1 || samplingTimeColumn == -1 || sampleCountColumn == -1 ){
            throw new RuntimeException("Not all required columns are present in the file");
        }

        String line = reader.readLine();

        while(line!=null){

            String[] entries = line.split(",");

            if (!idMap.containsKey("ID_"+entries[unitColumn])) {
                throw new RuntimeException("Trying to add a sampling event to unit " + entries[unitColumn] + " but" +
                        "this unit not previously defined");
            }

            InfectedUnit unit = idMap.get("ID_"+entries[unitColumn]);

            int sampCount = Integer.parseInt(entries[sampleCountColumn]);

            unit.addSamplingEvent(Double.parseDouble(entries[samplingTimeColumn]), sampCount);

            line = reader.readLine();
        }
    }

    // events are only relevant if there is a sampling event somewhere further up the tree

    private HashSet<FlexibleTree> makeTreelet(InfectedUnit unit, ArrayList<Event> relevantEvents, boolean forceCoalescence){

        if(relevantEvents.size()==0){
            return null;
        }

        ArrayList<SimpleNode> nodes = new ArrayList<SimpleNode>();

        unit.sortEvents();

        double lastRelevantEventTime = Double.NEGATIVE_INFINITY;

        for(Event event : relevantEvents){
            if(event.time > lastRelevantEventTime){
                lastRelevantEventTime = event.time;
            }
        }

        double activeTime = lastRelevantEventTime - unit.infectionEvent.time;

        for(Event event : relevantEvents){
            Taxon taxon;
            if(event.type == EventType.INFECTION){
                for(int instanceCount = 0; instanceCount < event.instances; instanceCount++){
                    taxon = new Taxon(event.infectee.id+"_infected_by_"+event.infector.id+"_"+event.time+"_lineage_"+(instanceCount+1));
                    taxon.setDate(new Date(event.time - unit.infectionEvent.time, Units.Type.YEARS, false));
                    SimpleNode node = new SimpleNode();
                    node.setTaxon(taxon);
                    nodes.add(node);
                    node.setHeight(unit.infectionEvent.time - event.time);
                    node.setAttribute("Event", event);
                }

            } else {
                for (int instanceCount = 0; instanceCount < event.instances; instanceCount++) {
                    taxon = new Taxon(unit.id + "_sampled_" + (instanceCount + 1) + "_" + event.time);
                    taxon.setDate(new Date(event.time - unit.infectionEvent.time, Units.Type.YEARS, false));
                    SimpleNode node = new SimpleNode();
                    node.setTaxon(taxon);
                    nodes.add(node);
                    node.setHeight(unit.infectionEvent.time - event.time);
                    node.setAttribute("Event", event);
                }
            }

        }

        HashSet<FlexibleNode> treeletRoots;

        if(nodes.size()>1){
            treeletRoots = simulateCoalescent(nodes, demFunct, activeTime, forceCoalescence);
        } else {
            treeletRoots = new HashSet<FlexibleNode>();
            FlexibleNode treeletRoot = new FlexibleNode(new SimpleTree(nodes.get(0)), nodes.get(0), true);
            treeletRoot.setHeight(0);

            treeletRoots.add(treeletRoot);
        }

        HashSet<FlexibleTree> out = new HashSet<FlexibleTree>();

        // add the root branch lengths

        for(FlexibleNode treeletRoot : treeletRoots){
            FlexibleNode infectionNode = new FlexibleNode();
            infectionNode.setHeight(activeTime);
            infectionNode.addChild(treeletRoot);
            treeletRoot.setLength(activeTime - treeletRoot.getHeight());
            infectionNode.setAttribute("Event", unit.infectionEvent);

            FlexibleTree aTree = new FlexibleTree(infectionNode);

            for(int i=0; i<aTree.getNodeCount(); i++){
                FlexibleNode node = (FlexibleNode)aTree.getNode(i);
                node.setAttribute("Unit", unit.id);
            }

            out.add(aTree);
        }



        return out;
    }

    private ArrayList<FlexibleTree> makeTrees(boolean forceCoalescences){

        // find the first case

        ArrayList<InfectedUnit> introducedCases = new ArrayList<InfectedUnit>();

        for(InfectedUnit unit : units){
            if(unit.parent==null){
                introducedCases.add(unit);
            }
        }

        if(introducedCases.size()==0){
            throw new RuntimeException("Can't find a first case");
        }

        ArrayList<FlexibleTree> out = new ArrayList<FlexibleTree>();

        for(InfectedUnit introduction : introducedCases) {
            if(introduction.childEvents.size()>0) {

                coalescentProbability = 1;

                System.out.println("Building tree(s) for descendants of " + introduction.id);
                ArrayList<FlexibleNode> outTreeRoots = makeSubtree(introduction, forceCoalescences);

                for(FlexibleNode outTreeRoot : outTreeRoots) {

                    FlexibleTree finalTree = new FlexibleTree(outTreeRoot, false, true);
                    finalTree.setAttribute("firstCase", introduction.id);
                    out.add(finalTree);

                }


                if(outTreeRoots.size()==0){
                    progressStream.println("This individual has no sampled descendants");
                }
                if (coalescentProbability < 0.9 & forceCoalescences) {
                    progressStream.println("WARNING: any phylogeny for descendants of " + introduction.id + " is quite " +
                            "improbable (p<" + (coalescentProbability) + ") given this demographic function. Consider " +
                            "another.");
                }

                System.out.println();
            }


        }
        return out;
    }

    // make the tree from this unit up

    private ArrayList<FlexibleNode> makeSubtree(InfectedUnit unit, boolean forceCoalescence){

        HashMap<Event, ArrayList<FlexibleNode>> eventToSubtreeRoots = new HashMap<>();

        ArrayList<Event> relevantEvents = new ArrayList<Event>();

        for(Event event : unit.childEvents){

            if(event.type == EventType.INFECTION){

                ArrayList<FlexibleNode> childSubtreeRoots = makeSubtree(event.infectee, forceCoalescence);

                event.setInstanceCount(childSubtreeRoots.size());

                if(childSubtreeRoots.size() != 0){
                    relevantEvents.add(event);
                    eventToSubtreeRoots.put(event, childSubtreeRoots);
                }

            } else if(event.type == EventType.SAMPLE) {
                relevantEvents.add(event);
            }
        }

        HashSet<FlexibleTree> unitTreelets = makeTreelet(unit, relevantEvents, forceCoalescence);

        if(unitTreelets.size()==0){
            return null;
        }


        for(Event event : relevantEvents){
            if (event.type == EventType.INFECTION) {
                ArrayList<FlexibleNode> relevantTips = new ArrayList<>();
                for (FlexibleTree unitTreelet : unitTreelets) {
                    for (int i = 0; i < unitTreelet.getExternalNodeCount(); i++) {
                        FlexibleNode tip = (FlexibleNode) unitTreelet.getExternalNode(i);

                        Event tipEvent = (Event) unitTreelet.getNodeAttribute(tip, "Event");

                        if (tipEvent == event) {
                            relevantTips.add(tip);
                        }
                    }
                }

                if (relevantTips.size() != eventToSubtreeRoots.get(event).size()) {
                    throw new RuntimeException("Numbers do not match");
                }


                for (int i = 0; i < relevantTips.size(); i++) {

                    FlexibleNode tip = relevantTips.get(i);
                    FlexibleNode root = eventToSubtreeRoots.get(event).get(i);

                    FlexibleNode firstSubtreeSplit = root.getChild(0);

                    root.removeChild(firstSubtreeSplit);
                    tip.addChild(firstSubtreeSplit);
                }
            }
        }

        ArrayList<FlexibleNode> out = new ArrayList<>();

        for(FlexibleTree unitTreelet : unitTreelets){
            out.add((FlexibleNode)unitTreelet.getRoot());
        }

        return out;
    }

    private HashSet<FlexibleNode> simulateCoalescent(ArrayList<SimpleNode> nodes, DemographicFunction demogFunct,
                                            double maxHeight, boolean forceCoalescence){

        double earliestNodeHeight = Double.NEGATIVE_INFINITY;

        for(SimpleNode node : nodes){
            if(node.getHeight()>earliestNodeHeight){
                earliestNodeHeight = node.getHeight();
            }
        }
        double maxLastInterval = earliestNodeHeight;
        double probNoCoalesenceInTime = Math.exp(demogFunct.getIntensity(maxLastInterval));

        coalescentProbability *= (1-probNoCoalesenceInTime);

        CoalescentSimulator simulator = new CoalescentSimulator();

        SimpleNode[] simResults;
        int failCount = 0;

        do {
            simResults = simulator.simulateCoalescent(nodes.toArray(new SimpleNode[nodes.size()]),
                    demogFunct, -maxHeight, 0, true);
            if(simResults.length>1){
                failCount++;
                System.out.println("Failed to coalesce lineages: "+failCount);
            }
        } while(simResults.length!=1 & forceCoalescence);

        HashSet<FlexibleNode> out = new HashSet<FlexibleNode>();

        for(int rootCount = 0; rootCount < simResults.length; rootCount++){
            SimpleNode root = simResults[rootCount];

            SimpleTree simpleTreelet = new SimpleTree(root);

            for (int i=0; i<simpleTreelet.getNodeCount(); i++) {
                SimpleNode node = (SimpleNode)simpleTreelet.getNode(i);
                node.setHeight(node.getHeight() + maxHeight);
            }

            out.add(new FlexibleNode(simpleTreelet, root, true));
        }

        return out;

    }

    private FlexibleTree makeWellBehavedTree(FlexibleTree tree){
        FlexibleTree newPhylogeneticTree = new FlexibleTree(tree, false);

        newPhylogeneticTree.beginTreeEdit();
        for(int i=0; i<newPhylogeneticTree.getInternalNodeCount(); i++){
            FlexibleNode node = (FlexibleNode)newPhylogeneticTree.getInternalNode(i);
            if(newPhylogeneticTree.getChildCount(node)==1){
                FlexibleNode parent = (FlexibleNode)newPhylogeneticTree.getParent(node);
                FlexibleNode child = (FlexibleNode)newPhylogeneticTree.getChild(node, 0);
                if(parent!=null){
                    double childHeight = newPhylogeneticTree.getNodeHeight(child);
                    newPhylogeneticTree.removeChild(parent, node);
                    newPhylogeneticTree.addChild(parent, child);
                    newPhylogeneticTree.setNodeHeight(child, childHeight);
                } else {
                    child.setParent(null);
                    newPhylogeneticTree.setRoot(child);
                }
            }
        }
        newPhylogeneticTree.endTreeEdit();


        return new FlexibleTree(newPhylogeneticTree, true);
    }

    private class InfectedUnit{
        private String id;
        private ArrayList<Event> childEvents;
        private Event infectionEvent;
        private InfectedUnit parent;

        private InfectedUnit(String id){
            this.id = id;
            parent = null;
            childEvents = new ArrayList<Event>();
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

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("virusTreeBuilder", "<infections-file-name> <sample-file-name> <output-file-name-root>");
    }



    public static void main(String[] args){

        ModelType model = ModelType.CONSTANT;
        double startNe = 1;
        double growthRate = 0;
        double t50 = 0;

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption(DEMOGRAPHIC_MODEL, demographics, false, "The type of within-host" +
                                " demographic function to use, default = constant"),
                        new Arguments.RealOption(STARTING_POPULATION_SIZE,"The effective population size at time zero" +
                                " (used in all models), default = 1"),
                        new Arguments.RealOption(GROWTH_RATE,"The effective population size growth rate (used in" +
                                " exponential and logistic models), default = 0"),
                        new Arguments.RealOption(T50,"The time point, relative to the time of infection in backwards " +
                                "time, at which the population is equal to half its final asymptotic value, in the " +
                                "logistic model default = 0"),
                        new Arguments.Option(FORCE_COALESCENCE,"Whether to insist on coalesence of all lineages before " +
                                "infection or allow incomplete bottlenecks"),
                        new Arguments.LongOption(SEED, "The random number seed")
                });


        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption(HELP)) {
            printUsage(arguments);
            System.exit(0);
        }

        if (arguments.hasOption(DEMOGRAPHIC_MODEL)) {
            String modelString = arguments.getStringOption(DEMOGRAPHIC_MODEL);
            if(modelString.toLowerCase().startsWith("c")){
                model = ModelType.CONSTANT;
            } else if(modelString.toLowerCase().startsWith("e")){
                model = ModelType.EXPONENTIAL;
            } else if(modelString.toLowerCase().startsWith("l")){
                model = ModelType.LOGISTIC;
            } else {
                progressStream.print("Unrecognised demographic model type");
                System.exit(1);
            }
        }


        if(arguments.hasOption(STARTING_POPULATION_SIZE)){
            startNe = arguments.getRealOption(STARTING_POPULATION_SIZE);
        }

        if(arguments.hasOption(GROWTH_RATE) && model!= ModelType.CONSTANT){
            growthRate = arguments.getRealOption(GROWTH_RATE);
        }

        if(arguments.hasOption(T50) && model== ModelType.LOGISTIC){
            t50 = arguments.getRealOption(T50);
        }

        if(arguments.hasOption(SEED)){
            MathUtils.setSeed(arguments.getLongOption(SEED));
        }

        boolean forceCoalesence = arguments.hasOption(FORCE_COALESCENCE);

        DemographicFunction demoFunction = null;

        switch(model){
            case CONSTANT: {
                demoFunction = new ConstantPopulation(Units.Type.YEARS);
                ((ConstantPopulation)demoFunction).setN0(startNe);
                break;
            }
            case EXPONENTIAL: {
                demoFunction = new ExponentialGrowth(Units.Type.YEARS);
                ((ExponentialGrowth)demoFunction).setN0(startNe);
                ((ExponentialGrowth)demoFunction).setGrowthRate(growthRate);
                break;
            }
            case LOGISTIC: {
                demoFunction = new LogisticGrowthN0(Units.Type.YEARS);
                ((LogisticGrowthN0)demoFunction).setN0(startNe);
                ((LogisticGrowthN0)demoFunction).setGrowthRate(growthRate);
                ((LogisticGrowthN0)demoFunction).setT50(t50);
                break;
            }
        }

        final String[] args2 = arguments.getLeftoverArguments();

        if(args2.length!=3){
            printUsage(arguments);
            System.exit(1);
        }

        String infectionsFileName = args2[0];
        String samplesFileName = args2[1];
        String outputFileRoot = args2[2];


        TransmissionTreeToVirusTree3 instance = new TransmissionTreeToVirusTree3(samplesFileName,
                infectionsFileName, demoFunction, outputFileRoot);

        try{
            instance.run(forceCoalesence);
        } catch (IOException e){
            e.printStackTrace();
        }
    }

}
