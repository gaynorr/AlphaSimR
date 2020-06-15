#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <queue>
//#include<stack>
#include <boost/weak_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/intrusive_ptr.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include "constants.h"

using namespace std;



// LIST OF TYPES DECLARED HERE

typedef unsigned int uint;

// sort by event time
struct byEventTime;
// sort by end position of a gene conversion tract
struct byEndPos;
// sort by allele freq in ascending order
struct byAlleleFreq;
// sort by the height of the bottom node of an edge
struct byBottomNode;
// sort by node pointer address
struct byNodePtr;

class ChrPosition;
typedef queue<ChrPosition> ChrPositionQueue;

class PtrRefCountable;
// A graph node
class Node;
// Smart pointer wrapper classes to automate garbage collection
//typedef boost::shared_ptr<Node> NodePtr;
typedef boost::intrusive_ptr<Node> NodePtr;
typedef vector <NodePtr> NodePtrVector;
typedef list <NodePtr> NodePtrList;
typedef set<NodePtr,byNodePtr> NodePtrSet;


class Edge;
typedef boost::shared_ptr<Edge> EdgePtr;
//typedef boost::intrusive_ptr<Edge> EdgePtr;
// Weak pointer that observes a shared pointer. Breaks cycles.
typedef boost::weak_ptr<Edge> WeakEdgePtr;
// For random access to any edge
typedef vector <EdgePtr> EdgePtrVector;
// Classify edges by population
typedef vector <EdgePtrVector> EdgePtrVectorByPop;
typedef list<EdgePtr> EdgePtrList;
// stores a collection of vector indices of elements that were deleted in the EdgePtrVector data structure
typedef queue<int> EdgeIndexQueue;
typedef vector <EdgeIndexQueue>  EdgeIndexQueueByPop;

class Event;
// shared pointer wrapper
typedef boost::intrusive_ptr<Event> EventPtr;
typedef list<EventPtr> EventPtrList;

class Population;
// vector of population objects used particularly
// in caching the pop counts in event traversal
typedef vector<Population> PopVector;

class GeneConversion;
typedef GeneConversion * GeneConversionPtr;
typedef set<GeneConversionPtr,byEndPos> GeneConversionPtrSet;

class HotSpotBin;
typedef HotSpotBin * HotSpotBinPtr;
typedef list<HotSpotBinPtr> HotSpotBinPtrList;

class AlleleFreqBin;
typedef AlleleFreqBin*  AlleleFreqBinPtr;
typedef set<AlleleFreqBinPtr,byAlleleFreq> AlleleFreqBinPtrSet;

class Mutation;
typedef Mutation * MutationPtr;
typedef vector<MutationPtr> MutationPtrVector;

typedef vector<vector<string> > CommandArguments;
typedef vector<vector<double> > MatrixDouble;


// TYPE DEFINITIONS BELOW

// This class can be called by boost's intrusive_ptr class to
// update the number of objects pointing to this object.
// Pointers of subclasses which are wrapped by intrusive_ptr
// should extend this class.  This should be less overhead
// than shared_ptr since the reference counter is not allocated
// from the heap again.


class PtrRefCountable
{
public:
  PtrRefCountable();
  // Declaring virtual ensures subclass destructor is called first.
  virtual ~PtrRefCountable();
protected:
  long references;
  // called by intrusive_ptr when to increment references
  friend void intrusive_ptr_add_ref(PtrRefCountable * p){
    ++(p->references);
  }
  // called by intrusive_ptr when to decrement references, deleting
  // object when necessary
  friend void intrusive_ptr_release(PtrRefCountable * p){
    if (--(p->references)==0) {
      delete p;
    }
  }
};

class ChrPosition
{
public:
  ChrPosition(unsigned long int iGraphIteration,double position);
  unsigned long int iGraphIteration;
  double position;
};

// The data structure encapsulating details about a population

class Population
{
public:
  Population();
  // the following accessor methods allows the object to perform sanity
  // checks on private variables (i.e. values must be positive)
  void setChrSampled(int);
  // get number of chromosomes remaining at a given time for this population
  int getChrSampled();
  // changes number of chromosomes by param value
  void changeChrSampled(int change);
  // assigns the population size relative to 1.
  void setPopSize(double dPopSize);
  // gets the population size
  double getPopSize();
  void setGrowthAlpha(double);
  double getGrowthAlpha();
  void setLastTime(double);
  double getLastTime();
private:
  double dLastTime;
  double dGrowthAlpha;
  int iChrSampled;
  double dPopSize;
};



// Edges connect graph nodes
class Edge:public PtrRefCountable
{
public:
  Edge(NodePtr & topNode,NodePtr & bottomNode);
  ~Edge();
  Edge();
  NodePtr & getTopNodeRef();
  NodePtr & getBottomNodeRef();
  double getLength();
  
  // If you call the following two methods, you better know what you're doing.
  // Any sets ordered by top or bottom node height
  // can be misordered.
  //    void setTopNode(NodePtr & topNode);
  void setBottomNode(NodePtr & bottomNode);
  // Mark an edge that is deleted since removing from vectors is costly
  bool bDeleted;
  // Mark an edge that is deleted that it is queue so it can be re-used
  bool bInQueue;
  // Denotes whether edge is part of the current local tree
  bool bInCurrentTree;
  // The current iteration along the unit length chromosome
  int iGraphIteration;
private:
  NodePtr topNode,bottomNode;
  double dLength;
};





class Node:public PtrRefCountable
{
public:
  // constants for identifying this node type
  enum NodeType {COAL,XOVER,MIGRATION,SAMPLE,QUERY};
  // constants used when identifying whether to
  // invoke operations on top edges or bottom edges
  enum EdgeLocation {TOP_EDGE,BOTTOM_EDGE};
  
  Node(NodeType iType,short int iPopulation,double dHeight);
  ~Node();
  // returns the population for this sampled node, any line
  // going back in time from this node belongs to this population
  short int getPopulation();
  // an integer constant identifying this node. see constants below
  NodeType getType();
  // an integer constant identifying this node. see constants below
  const char * getTypeStr();
  // height in graph relative to time 0: present time
  double getHeight();
  // clear top edges
  void clearTopEdges();
  // clears bottom edges;
  //    void clearBottomEdges();
  short unsigned int getTopEdgeSize();
  short unsigned int getBottomEdgeSize();
  // a vector of the edge(s) above the node
  EdgePtr  getTopEdgeByIndex(short unsigned int index);
  EdgePtr  getBottomEdgeByIndex(short unsigned int index);
  // method to be called when connecting a new edge to this node
  // when an old edge is known to be replaced
  void replaceOldWithNewEdge(EdgeLocation iLocation,
                             EdgePtr & oldEdge,EdgePtr & newEdge);
  // method to be called when connecting a new edge to this node
  void addNewEdge(EdgeLocation iLocation,EdgePtr & newEdge);
  // invalidates the event associated with this node.
  void removeEvent();
  // is there an event associated with this node?
  bool isEventDefined() {
    return false;
  }
  
  // assign an event with this node
  void setEvent(EventPtr & assoEvent);
  // a place holder to identify a node at the top of the coalescing
  // line.  the height of this node assures that the coalescing
  // line will always be included as a candidate for coalescence
  // when traversing event list
  static const double MAX_HEIGHT;
  bool bDeleted;
  
private:
  EventPtr event;
  WeakEdgePtr topEdge1,topEdge2,bottomEdge1,bottomEdge2;
  //    EdgePtr topEdge1,topEdge2,bottomEdge1,bottomEdge2;
  
  bool bEventDefined;
  short int iPopulation;
  NodeType iType;
  short unsigned int topEdgeSize,bottomEdgeSize;
  double dHeight;
};

//double Node::MAX_HEIGHT;



class SampleNode:public Node
{
public:
  SampleNode(short int iPopulation,int iId);
  int iId;
  bool bAffected;
};


class Event:public PtrRefCountable
{
public:
  enum EventType {
    PAST_COAL,
    NEW_COAL,
    GROWTH,
    MIGRATION_RATE,
    MIGRATION_MATRIX_RATE,
    XOVER,
    NEW_MIGRATION,
    PAST_MIGRATION,
    GLOBAL_POPSIZE,
    GLOBAL_POPGROWTH,
    GLOBAL_MIGRATIONRATE,
    POPSIZE,
    POPSPLIT,
    POPJOIN};
  // everyt event must have a type and time associated with it
  Event(EventType iType,double dTime);
  ~Event();
  // returns the type based on the constants below
  EventType getType();
  //    int iGraphIteration;
  // time on the graph from time 0 at bottom
  double getTime();
  bool bMarkedForDelete;
  
private:
  EventType iType;
  double dTime;
};

// For simple events.

class GenericEvent:public Event
{
public:
  GenericEvent(EventType iType,double dTime,double dParameterValue);
  double getParamValue();
private:
  double dParameterValue;
};

// Represents events such as change in pop size or a new pop growth

class PopSizeChangeEvent:public Event
{
public:
  PopSizeChangeEvent(EventType iType,double dTime,short int iPopulationIndex,
                     double dPopChangeParam);
  short int getPopulationIndex();
  double getPopChangeParam();
private:
  short int iPopulationIndex;
  double dPopChangeParam;
};

// Represents actual migration events (one chr from one pop to another pop)

class MigrationEvent:public Event
{
public:
  MigrationEvent(EventType iType,double dTime,
                 short int iPopMigratedTo,short int iPopMigratedFrom);
  short int getPopMigratedTo();
  short int getPopMigratedFrom();
private:
  short int iPopMigratedTo,iPopMigratedFrom;
};

// Two populations merge

class PopJoinEvent:public Event
{
public:
  PopJoinEvent(EventType iType,double dTime,
               short int iSourcePop,short int iDestPop);
  short int getSourcePop();
  short int getDestPop();
private:
  short int iSourcePop,iDestPop;
};

// A change in migration rate

class MigrationRateEvent:public Event
{
public:
  MigrationRateEvent(EventType iType,double dTime,
                     short int iSourcePop,short int iDestPop,
                     double dRate);
  short int getSourcePop();
  short int getDestPop();
  double getRate();
private:
  short int iSourcePop,iDestPop;
  double dRate;
};



// A new configuration where the entire matrix is used

class MigrationRateMatrixEvent:public Event{
public:
  MigrationRateMatrixEvent(EventType iType,double dTime,
                           MatrixDouble dMigrationMatrix);
  MatrixDouble getMigrationMatrix();
  //    ~MigrationRateMatrixEvent();
private:
  MatrixDouble dMigrationMatrix;
};

// A coalescent event

class CoalEvent:public Event
{
public:
  CoalEvent(EventType iType,double dTime,
            short int iPopulation);
  short int getPopulation();
  
private:
  short int iPopulation;
};

// A crossover event

class XoverEvent:public Event
{
public:
  XoverEvent(EventType iType,double dTime,
             short int iPopulation);
  short int getPopulation();
private:
  short int iPopulation;
};


class GeneConversion:public PtrRefCountable
{
public:
  GeneConversion(double dEndPos);
  ~GeneConversion();
  NodePtr xOverNode;
  double dEndPos;
};

//typedef boost::intrusive_ptr<GeneConversion> GeneConversionPtr;






// a bin for variable recombination rates

class HotSpotBin{
public:
  HotSpotBin(double dStart,
             double dEnd,
             double dRate);
  double dStart;
  double dEnd;
  double dRate;
};

//typedef boost::shared_ptr<HotSpotBin> HotSpotBinPtr;


// a bin for an allele frequency ascertainment scheme

class AlleleFreqBin{
public:
  AlleleFreqBin(double dStart,
                double dEnd,double dFreq);
  ~AlleleFreqBin();
  
  double dStart;
  double dEnd;
  double dFreq;
  int iObservedCounts;
  
};

class AlphaSimRReturn {
public:
  AlphaSimRReturn();
  vector<bool > haplotypes;
  double length;
};

class Mutation{
public:
  Mutation(double dLocation,double dFreq);
  ~Mutation();
  double dFreq;
  double dLocation;
  bool bPrintOutput;
};

// Configuration container populated by parameter reading procedure
// can be used by any simulator implementation
class Configuration
{
public:
  double dTheta,dGlobalMigration,dRecombRateRAcrossSites,
  dGeneConvRatio,dSeqLength,dBasesToTrack;
  unsigned int iSampleSize,iIterations,iGeneConvTract;
  unsigned short int iTotalPops;
  long iRandomSeed;
  bool bSNPAscertainment,bFlipAlleles;
  //    bool bHighMutationRate;
  bool bVariableRecomb,bNewickFormat;
  bool bMigrationChangeEventDefined;
  MatrixDouble dMigrationMatrix;
  // stores population info
  PopVector pPopList;
  // stores the user defined events
  EventPtrList * pEventList;
  HotSpotBinPtrList * pHotSpotBinPtrList;
  AlleleFreqBinPtrSet * pAlleleFreqBinPtrSet;
  Configuration();
  ~Configuration();
};

class RandNumGenerator{
public:
  RandNumGenerator(unsigned long iRandomSeed);
  ~RandNumGenerator();
  // Given the parameter lambda, returns an exponentially distributed
  // random variable
  double expRV(double dLambda);
  // Returns a uniforumly distributed random variable
  double unifRV();
private:
  boost::uniform_01<boost::mt19937> * unif;
};


// The main class that implements the Wiuf and Hein
// algorithm, incorporating demographic events.


// Graph based implementation of a coalescent simulator.
class GraphBuilder
{
public:
  // This object is initialized with configuration supplied by the user
  GraphBuilder(Configuration *,RandNumGenerator *);
  ~GraphBuilder();
  // The entry point for building the graph while traversing the
  // the chromosome on the unit interval.
  void build();
  // Print the haplotypes in MS format
  void printHaplotypes();
  vector<AlphaSimRReturn> getMutations();
  
private:
  // The random number generator
  RandNumGenerator *pRandNumGenerator;
  // Points to user specified parameters
  Configuration *pConfig;
  
  vector<AlphaSimRReturn> mutations;
  // *** ESSENTIAL CONTAINERS POINTING TO
  // EDGES IN THE GRAPH AND THE MRCAS
  // a linked list of edges on the ARG
  EdgePtrList *pEdgeListInARG;
  // the vector containing all the edges of the last local tree
  EdgePtrVector * pEdgeVectorInTree;
  // the total edges that are valid in pEdgeVectorInTree
  // we update and access this value instead of calling size()
  // on pEdgeVectorInTree so that we avoid a costly allocation
  // and de allocation of the local tree list at every graph
  // iteration.
  unsigned int iTotalTreeEdges;
  // total branch length of the ARG
  double dArgLength;
  // Contains the total branch length of all edges
  // Should be updated when necessary by initializeCurrentTree()
  double dLastTreeLength;
  // grandMRCA stores the Node object which is the MRCA of the ARG
  // localMRCA is used when searching for the MRCA of the last tree
  NodePtr grandMRCA,localMRCA,xOverNode;
  constexpr static const auto dEpsilon = 1e-6;
  
  // *** USED FOR COMPUTING TIME TO COALESCENCE
  // linked list of events that must be traversed in order
  // to get accurate snapshots of the demographic parameters
  // at a given time.
  EventPtrList * pEventList;
  // matrix of migration rates
  MatrixDouble dMigrationMatrix;
  // the running edge that runs upwards from the xover point
  // used by the event traversing method to insert migration events
  EdgePtr coalescingEdge;
  // the running edge that runs upwards from the previous MRCA
  // used by the event traversing method to insert migration events
  EdgePtr originExtension;
  
  
  // *** USED FOR SEEKING EDGES TO COALESCE
  // A vector by pop of vector of the edges in the graph
  // for random access when picking a random spot
  // for coalescing to
  EdgePtrVectorByPop *pEdgeVectorByPop;
  // A vector by pop of a queue of integers to contain
  // the index in pEdgeVector that has been freed up from
  // a past deletion
  EdgeIndexQueueByPop *pVectorIndicesToRecycle;
  
  // Store the mutations so site ascertainment can be carried out in RAM
  MutationPtrVector * pMutationPtrVector;
  // workspace: array of affection status used for computing allele frequencies
  bool *  sites;
  
  // an array that stores the edge that have
  // yet to find the MRCA when constructing the last tree
  EdgePtr * pTreeEdgesToCoalesceArray;
  // vector containing the samples at time 0
  NodePtr * pSampleNodeArray;
  
  // recombination rate scaled to account for gene conversion to xover ratio.
  double dScaledRecombRate;
  // the iteration ID of the algorithm, corresponding to the
  // number of recombination events as spelled out by Wiuf and Hein algorithm
  int iGraphIteration;
  // stores the positions on the chromosome
  ChrPositionQueue * pChrPositionQueue;
  double dTrailingGap;
  bool bIncrementHistory;
  
  
  // *** THE FOLLOWING ARE USED FOR GENE CONVERSION
  // a set of gene conversion objects ordered by the end point of the tract.
  GeneConversionPtrSet * pGeneConversionPtrSet;
  // did we just start a gene conversion event?
  bool bBeginGeneConversion;
  // flag to close a pending gene conversion event
  bool bEndGeneConversion;
  // if gene conversion is to be closed at this iteration use the following
  // two saved edges
  EdgePtr gcOldEdge,gcNewEdge;
  
  
  // The following are the primary methods used in the build() method.
  
  // this function, loosely based on Hudson's algorithm, serves two purposes
  // 1. Initializes the graph(tree) from present time (0) to the grand ancestor.  Nodes
  // (coalescent and migration) and edges are generated; coalescent and
  // migration events are added to the existing event list specified by the
  // user.  A vector of Population objects are used to track changes
  // in sample size and pop size as time moves back.  A migration matrix
  // is also used to track changes in migration rates between populations
  // as we move up the tree.
  // 2. The graph is modified after a recombination event. A new time for
  // coalescence is computed by traversing the event list, adding new
  // events/nodes/edges along the way if necessary. When a new line extending
  // from the recombination node is eligible for coalescence, the number of
  // edges at the current time matching the recombination nodes's population is
  // considered as the rate and used as intensity parameter for an
  // exponential draw to select a proposed coalescent time.
  void traverseEvents(bool bBuildFromEventList,
                      NodePtr & xOverNode, EventPtr & newCoalEvent);
  // Modifies the graph after finding the next position on the chr where
  // a xover occurs.  Here a uniform draw on the total branch length
  // selects the branch and the location the xover occus on upon the graph.
  // To compute time to coalescence for the new lineage above the recombination
  // node, the function relies on traverseEvents()
  void invokeRecombination(GeneConversionPtr &);
  // check for any pending gene conversions to close at current position
  // If we find that the end point of an existing gene conversion is less
  // the current position, we need to backtrack to that end point.
  bool checkPendingGeneConversions(double & curPos);
  // Add a mutation uniformly to the local tree, allowing it to trickle down
  // to the sampled chromosomes
  void addMutations(double startPos,double endPos);
  // Uniform randomly selects an edge (and position) to insert a xover or mutation node into
  EdgePtr getRandomEdgeOnTree(double & dSplitPoint,double dRandSpot);
  // Once a coalescent height is determined, uniform randomly select an
  // eligible existing edge to coalesce with
  EdgePtr getRandomEdgeToCoalesce(EdgePtr & coalescingEdge,
                                  double dHeight);
  // Mark the edges within the ARG that corresponds to the local tree
  // at the current position
  void markCurrentTree();
  // Prune all edges with tree histories less than threshold specified by user
  void pruneARG(int iHistoryMax);
  // get the next pos based on hot spot rates
  bool getNextPos(double & curPos,HotSpotBinPtrList::iterator & hotSpotIt);
  // get the rate until next xover or gene conversion
  double getRate();
  
  // The following methods are utility methods called by the primary methods
  // above.
  
  // Retallies the total edge length of the ARG.  Is called after a
  // series of branches are added or deleted
  //    void InitializeGraphLength();
  // Retallies the total edge length of the last tree.  Is called
  // after a series of branches are added or deleted
  void initializeCurrentTree();
  // Traverses recusively down all edges until a sample node is reached
  // where the mutation bit is set.
  void mutateBelowEdge(EdgePtr & edge);
  // color edges above. bFirstSample is true only
  // the first time this method is called so all following calls
  // can try to find an grandMRCA at or greater than that found
  // from the first call. bCalledFromParent is true only
  // when called from MarkCurrentTree and false when
  // called recursively
  bool markEdgesAbove(bool bFirstSample,bool bCalledFromParent,
                      EdgePtr & edgeAboveNode,unsigned int & iSampleIndex);
  // recursively traverses graph until a coal or sample node
  // is reached.  This is called when a crossover is invoked,
  // so that the descendants of the crossover node are notified that
  // they are part of the most current tree.
  void markEdgesBelow(EdgePtr & edgeBelowNode);
  // Traverses recursively up migration nodes until a coalescent event
  // is reached in which the two other joining edges are joined
  void pruneEdgesAbove(EdgePtr & edgeAboveNode);
  // Traverses recusively down migration nodes until a coalescent event
  // is reached. This is used when an old higher MRCA is to be removed.
  void pruneEdgesBelow(EdgePtr & edgeBelowNode);
  
  
  // Mark edge as deleted and delete edge from all data structures
  // when efficient
  void deleteEdge(EdgePtr & edge);
  // Insert edge into all data structures pertaining to ARG
  void addEdge(EdgePtr & edge);
  // Add to current tree
  void addEdgeToCurrentTree(EdgePtr & edge);
  // Removes a node from an edge
  void mergeEdges(EdgePtr & topEdge,EdgePtr & bottomEdge);
  // Given a single edge, splits the edge and places a node in between
  void insertNodeInEdge(NodePtr & newNode,EdgePtr & oldEdge);
  // Given a temporary edge (the coalescing or grandMRCA extension edge),
  // splits the edge and places a node in between
  void insertNodeInRunningEdge(NodePtr & newNode,EdgePtr & tempEdge);
  // Dumps the current data structures to stdout for debugging
  void printDataStructures();
  // Checks that at a given time, the running totals for each population
  // match the edge counts at the level in the graph
  // Expensive and only called in debugging
  void checkPopCountIntegrity(PopVector &,double dTime);
  // Prints local tree in Newick format as done in MS
  string getNewickTree(double lastCoalHeight,NodePtr & curNode);
  
  
};

// The entry point for the program.
class Simulator
{
public:
  //*******************************************************
  // Creates a Configuration object by reading parameters
  // from the command line
  // We can always do more error checking here!
  void readInputParameters(CommandArguments args);
  // Calls any coalescent simulator (e.g. fastcoal, MS). In this
  // case, constructs a new graphbuilder and calls the build() function
  void beginSimulation();
  vector<AlphaSimRReturn> beginSimulationMemory();
  Simulator();
  ~Simulator(); //destructor
  
private:
  // Contains configuration information provided by the user.
  Configuration *pConfig;
};

// IMPLEMENTATIONS ARE BELOW

// SPECIAL SORT OPERATORS FOR STL CONTAINERS ARE IMPLEMENTED HERE

struct byNodePtr{
  bool operator()(const NodePtr& node1, const NodePtr& node2) const{
    return (node1<node2);
  }
};


struct byEventTime{
  bool operator()(const EventPtr& event1, const EventPtr& event2) const{
    return (event1->getTime()<event2->getTime());
  }
};

// sort by bottom node of an edge
struct byBottomNode{
  bool operator()(const EdgePtr& edge1, const EdgePtr& edge2) const{
    return (edge1->getBottomNodeRef()->getHeight()>edge2->getBottomNodeRef()->getHeight());
  }
};

struct byEndPos{
  bool operator()(const GeneConversionPtr& gc1, const GeneConversionPtr& gc2) const{
    return (gc1->dEndPos<gc2->dEndPos);
  }
};

struct byAlleleFreq{
  bool operator()(const AlleleFreqBinPtr & bin1, const AlleleFreqBinPtr & bin2) const{
    return (bin1->dStart<bin2->dStart && bin1->dEnd<bin2->dEnd);
  }
};


// INLINE FUNCTIONS ARE IMPLEMENTED HERE


inline int Population::getChrSampled(){
  return this->iChrSampled;
}

inline void Population::changeChrSampled(int change){
  //setChrSampled(this->iChrSampled + change);
  this->iChrSampled+=change;
}

inline double Population::getPopSize(){
  return this->dPopSize;
}

inline NodePtr & Edge::getTopNodeRef(){
  return this->topNode;
}

inline NodePtr & Edge::getBottomNodeRef(){
  return this->bottomNode;
}

inline double Edge::getLength(){
  return this->dLength;
}


inline short int Node::getPopulation(){
  return this->iPopulation;
}

// an enumeration identifying this node. see constants below
inline Node::NodeType Node::getType(){
  return this->iType;
}
// height in graph relative to time 0: present time
inline double Node::getHeight(){
  return this->dHeight;
}
// clear top edges
inline void Node::clearTopEdges(){
  this->topEdgeSize = 0;
}

inline short unsigned int Node::getTopEdgeSize(){
  return this->topEdgeSize;
}

inline short unsigned int Node::getBottomEdgeSize(){
  return this->bottomEdgeSize;
}

// a vector of the edge(s) above the node
inline EdgePtr Node::getTopEdgeByIndex(short unsigned int index){
  return index != 0u ?topEdge2.lock():topEdge1.lock();
  //        return index?topEdge2:topEdge1;
  //    return edge;
}

inline EdgePtr Node::getBottomEdgeByIndex(short unsigned int index){
  return index?bottomEdge2.lock():bottomEdge1.lock();
  //        return index?bottomEdge2:bottomEdge1;
  //    return edge;
}

inline void Node::removeEvent(){
  if (this->bEventDefined){
    this->event->bMarkedForDelete = true;
  }
}

inline void Node::setEvent(EventPtr & event){
  this->bEventDefined = true;
  this->event = event;
}

inline double Event::getTime(){
  return this->dTime;
}

inline Event::EventType Event::getType(){
  return this->iType;
}

inline short int MigrationEvent::getPopMigratedTo(){
  return iPopMigratedTo;
}

inline short int MigrationEvent::getPopMigratedFrom(){
  return iPopMigratedFrom;
}

inline short int CoalEvent::getPopulation(){
  return iPopulation;
}

inline short int XoverEvent::getPopulation(){
  return iPopulation;
}

inline double RandNumGenerator::expRV(double dLambda){ // generates exponential random variables
  return -log(1.-(*unif)())/(dLambda);
}

inline double RandNumGenerator::unifRV() {
  return (*unif)();
}

inline double GraphBuilder::getRate(){
  return dLastTreeLength*dScaledRecombRate;
}
