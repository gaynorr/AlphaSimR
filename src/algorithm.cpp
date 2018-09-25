#include <Rcpp.h>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <time.h>
#include "simulator.h"

using namespace std;


EdgePtr GraphBuilder::getRandomEdgeOnTree(double & dSplitPoint,
                                          double dRandomSpot){
  if (this->bEndGeneConversion && dSplitPoint!=-1.){
    dSplitPoint = gcNewEdge->getBottomNodeRef()->getHeight()+0.;
    return this->gcNewEdge;
  }
  bool found = false;
  double dRunningLength = 0.0;
  EdgePtrVector::iterator it = pEdgeVectorInTree->begin();
  EdgePtr curEdge;
  // the linked list of tree edges are lined up so that a uniform point
  // on the tree can be selected for where the crossover occurs
  unsigned int counter=0;
  while(!found && counter<iTotalTreeEdges){
    curEdge = *it;
    if (!curEdge->bDeleted){
      if (dRandomSpot<dRunningLength+curEdge->getLength()){
        // the random spot is within the current segment
        // the reference var splitpoint returns the position
        // relative to the bottom of this edge that the xover
        // occurs
        dSplitPoint = curEdge->getBottomNodeRef()->getHeight()+
          (dRandomSpot - dRunningLength);
        found = true;
      }else{
        dRunningLength+=curEdge->getLength();
        //if (pConfig->bDebug){
        //Rcpp::Rcerr<<"getRandomEdge: running length: "<<dRunningLength<<", random spot "<<dRandomSpot<<"\n"; 
        //}
      }
    }else{
      //if (pConfig->bDebug){
      //Rcpp::Rcerr<<"getRandomEdge: Sorry edge was deleted\n"; 
      //}
    }
    ++counter;
    ++it;
  }
  if (!found) throw "RandomSpot was out of range for xover";
  return curEdge;
}

EdgePtr GraphBuilder::getRandomEdgeToCoalesce(EdgePtr & coalescingEdge,
                                              double dHeight){
  
  if (this->bEndGeneConversion){
    return this->gcOldEdge;
  }
  int iPopulation = coalescingEdge->getBottomNodeRef()->getPopulation();
  EdgePtrVector & pEdgeVector = this->pEdgeVectorByPop->at(iPopulation);
  
  EdgeIndexQueue & pVectorIndicesToRecycle = this->pVectorIndicesToRecycle->at(iPopulation);
  unsigned int edgeVectorSize = pEdgeVector.size();
  if (pConfig->bDebug){
    Rcpp::Rcerr<<"getRandomEdgeToCoalesce: selecting population "<<iPopulation<<" with avail edges: "<< edgeVectorSize<<"\n";
  }
  bool found = false;
  int iUnifPick =-1;
  int out_of_range = 0,already_deleted=0;
  while(!found){
    // guessing by randomly selecting elements from the vector of graph
    // edges appears to work pretty well here
    iUnifPick = static_cast<int>(pRandNumGenerator->unifRV()*
      edgeVectorSize);
    EdgePtr edge = pEdgeVector[iUnifPick];
    if (edge->bDeleted){
      // if this is marked for deletion, take the opportunity
      // to store the vector index in a queue so this element can
      // recycled when a new edge is to be added
      if (!edge->bInQueue){
        pVectorIndicesToRecycle.push(iUnifPick);
        edge->bInQueue = true;
      }
      if (pConfig->bDebug){
        //Rcpp::Rcerr<<"getRandomEdgeToCoalesce: edge already deleted \n";
        ++already_deleted;
      }
      // if the proposed coalescing time is within the endpoints of this
      // randomly picked edge, select this edge
    }else if (edge->getBottomNodeRef()->getHeight()<dHeight
                && edge->getTopNodeRef()->getHeight()>dHeight) {
      found = true;
    }else{
      if (pConfig->bDebug){
        //Rcpp::Rcerr<<"getRandomEdgeToCoalesce: Edge out of range\n";
        ++out_of_range;
      }
    }
  }
  EdgePtr & edge = pEdgeVector[iUnifPick];
  if (pConfig->bDebug){
    Rcpp::Rcerr<<"Stats: out of range edges: "<<out_of_range<<" and already deleted: "<<already_deleted<<endl;
    Rcpp::Rcerr<<"Edge to coalesce: "<<edge->getBottomNodeRef()->getHeight()<<" to "<<
      edge->getTopNodeRef()->getHeight()<<endl;
  }
  
  return edge;
}

void GraphBuilder::mutateBelowEdge(EdgePtr & edge){
  NodePtr & bottomNode = edge->getBottomNodeRef();
  if (bottomNode->getType()==Node::SAMPLE  ){
    SampleNode* sample = static_cast<SampleNode*>(bottomNode.get());
    // mark this sampled node as affected.
    sample->bAffected = true;
  }else{
    for (unsigned int i=0;i<bottomNode->getBottomEdgeSize();++i){
#ifdef DIAG
      if (bottomNode->getBottomEdgeByIndex(i)->bDeleted){
        throw "Mutate below edge: this edge should have been removed";
      }
#endif
      EdgePtr edge = bottomNode->getBottomEdgeByIndex(i);
      if (edge->iGraphIteration==iGraphIteration){
        mutateBelowEdge(edge);
      }
    }
  }
}

bool GraphBuilder::markEdgesAbove(bool bFirstSample,bool  bCalledFromParent,
                                  EdgePtr & selectedEdge, unsigned int & iSampleIndex){
  if (bFirstSample || bCalledFromParent || !selectedEdge->bInCurrentTree){
    if (!selectedEdge->bInCurrentTree) addEdgeToCurrentTree(selectedEdge);
    NodePtr topNode = selectedEdge->getTopNodeRef();
    double dTopHeight = topNode->getHeight();
    double dproposedOriginHt = localMRCA->getHeight();
    if (dTopHeight>=dproposedOriginHt &&
        topNode->getType()==Node::COAL){
      // If we reach at least as high as proposed MRCA
      // save this edge as the vector
      // so we don't traverse from the bottom
      pTreeEdgesToCoalesceArray[iSampleIndex]=selectedEdge;
      if (bFirstSample){
        if (dTopHeight>dproposedOriginHt)
          // establish a new proposed MRCA
          // for susequent samples to try reaching
          localMRCA = topNode;
        // proceed to the next sampled node
        return true;
      }else{
        if (dTopHeight>dproposedOriginHt){
          // if we surpassed a MRCA height previously proposed
          // upgrade the proposed MRCA height
          localMRCA = topNode;
          // start over in the loop to sample 0
          return false;
        }else{
          if (localMRCA!=topNode){
            Rcpp::Rcerr<<"proposed grandMRCA != top Node\n";
          }
          // here the current height is EQUAL to the proposed MRCA,
          // this is good and we can proceed to the next sampled node
          return true;
        }
      }
    }
    else{
      // if we come across a crossover, there are two edges above,
      // we choose the one with the more recent (larger) tree label
      if (topNode->getType()==Node::XOVER){
        EdgePtr edge0 =topNode->getTopEdgeByIndex(0);
        EdgePtr edge1 =topNode->getTopEdgeByIndex(1);
        int iter0=-1,iter1=-1;
        if (edge0!=NULL) iter0=edge0->iGraphIteration;
        if (edge1!=NULL) iter1=edge1->iGraphIteration;
        if (iter0==-1 && iter1==-1) throw "In mark edges above, both edges above xover were missing";
        EdgePtr & nextEdge = (iter0>iter1)? edge0:edge1;
        return markEdgesAbove(bFirstSample,false,nextEdge,
                              iSampleIndex);
        // for other nodes, just choose the single edge above it
        // recursively call this method again
      }else{
        EdgePtr edge = topNode->getTopEdgeByIndex(0);
        return markEdgesAbove(bFirstSample,false,edge,
                              iSampleIndex);
      }
    }
  }else{
    // the selected edge was in the current tree already and this is for
    // a sampled node after the first one and was called by this method.
    // the latter condition implies that we previously explored up this
    // node.  in this case, save time and proceed to the next node
    return true;
  }
}

void GraphBuilder::markEdgesBelow(EdgePtr & selectedEdge){
  selectedEdge->iGraphIteration = iGraphIteration;
  NodePtr  & bottomNode = selectedEdge->getBottomNodeRef();
  if (bottomNode->getType()==Node::COAL||
      bottomNode->getType()==Node::SAMPLE) return;
  else{
    EdgePtr bottomEdge = bottomNode->getBottomEdgeByIndex(0);
    markEdgesBelow(bottomEdge);
  }
  
}

void GraphBuilder::pruneEdgesBelow(EdgePtr & selectedEdge){
  NodePtr & bottomNode = selectedEdge->getBottomNodeRef();
  if (bottomNode->getType()==Node::COAL){
    // this is the terminating condition and should be the new
    // lower MRCA
    grandMRCA = bottomNode;
    grandMRCA->clearTopEdges();
  }else{
    // we didn't find a coalescent node (MRCA) so keep searching down
    bottomNode->removeEvent();
    EdgePtr bottomEdge = bottomNode->getBottomEdgeByIndex(0);
    pruneEdgesBelow(bottomEdge);
  }
  deleteEdge(selectedEdge);
}

void GraphBuilder::pruneEdgesAbove(EdgePtr & selectedEdge){
  NodePtr & topNode = selectedEdge->getTopNodeRef();
  short int topNodeType = topNode->getType();
  if (topNodeType==Node::COAL){
    EdgePtr edge0 = topNode->getBottomEdgeByIndex(0);
    EdgePtr edge1 = topNode->getBottomEdgeByIndex(1);
    if (topNode->getTopEdgeSize()==0){
#ifdef DIAG
      if (topNode!=grandMRCA) throw "Top node is NOT grandMRCA!\n";
#endif
      EdgePtr & otherEdge = edge0==selectedEdge ?
      edge1 : edge0;
      pruneEdgesBelow(otherEdge);
    }else{
      // just some coal node below the grandMRCA. fuse the other two edges.
      EdgePtr edgeAboveNode = topNode->getTopEdgeByIndex(0);
      EdgePtr & edgeBelowNode = edge0==
        selectedEdge?edge1:edge0;
      mergeEdges(edgeAboveNode,edgeBelowNode);
    }
    
  }else{
#ifdef DIAG
    if (topNodeType!=Node::MIGRATION
          && topNodeType!=Node::QUERY){
      Rcpp::Rcerr<<"At edge with "<<selectedEdge->getBottomNodeRef()->getHeight()<<
        " and "<<selectedEdge->getTopNodeRef()->getHeight()<<endl;
      throw "Prune up did not find a migration node";
    }
#endif
    EdgePtr topEdge = topNode->getTopEdgeByIndex(0);
    pruneEdgesAbove(topEdge);
  }
  topNode->removeEvent();
  deleteEdge(selectedEdge);
}

void GraphBuilder::invokeRecombination(GeneConversionPtr & geneConversionPtr){
  double dSplitPoint=0.0;
  double dRandomSpot = pRandNumGenerator->unifRV() * dLastTreeLength;
  EdgePtr selectedXoverEdge = getRandomEdgeOnTree(dSplitPoint,dRandomSpot);
  xOverNode = NodePtr(new Node(Node::XOVER,
                               selectedXoverEdge->getBottomNodeRef()->getPopulation(),dSplitPoint));
  // save old edge if gene conversion
  if (this->bBeginGeneConversion) geneConversionPtr->xOverNode = xOverNode;
  
  if(pConfig->bDebug){
    Rcpp::Rcerr<<"Adding xover point at "<<dSplitPoint<<
      " for edge "<<selectedXoverEdge->getBottomNodeRef()->getHeight()<<
        " to "<<selectedXoverEdge->getTopNodeRef()->getHeight()<<
          " of population "<<selectedXoverEdge->getBottomNodeRef()->getPopulation()<<endl;
  }
  EventPtr pXoverEventWrapper = EventPtr(new XoverEvent
                                           (Event::XOVER,dSplitPoint,xOverNode->getPopulation()));
  // insert the xover node at the selected edge
  insertNodeInEdge(xOverNode,selectedXoverEdge);
  
  // mark the edge above the xover.  this is the old coalescing line to prune.
  // set up a coalescing line that is visible to the
  // event building method
  NodePtr nodeAboveXover = NodePtr(new Node
                                     (Node::QUERY,xOverNode->getPopulation(),Node::MAX_HEIGHT));
  coalescingEdge = EdgePtr(new Edge(nodeAboveXover,xOverNode));
  nodeAboveXover->addNewEdge(Node::BOTTOM_EDGE,coalescingEdge);
  xOverNode->addNewEdge(Node::TOP_EDGE,coalescingEdge);
  // similarly, set up a line that runs straight up from the grandMRCA
  // in case migration events need to be inserted above the grandMRCA.
  NodePtr nodeAboveOrigin = NodePtr(new Node
                                      (Node::QUERY,grandMRCA->getPopulation(),Node::MAX_HEIGHT));
  originExtension = EdgePtr(new Edge(nodeAboveOrigin,grandMRCA));
  nodeAboveOrigin->addNewEdge(Node::BOTTOM_EDGE,originExtension);
  grandMRCA->addNewEdge(Node::TOP_EDGE,originExtension);
  // now initialize the parameters that will be sent to the event traversal routine
  EventPtr newCoalEvent;
  // Traverse the events to find the proper coalescent waiting time.
  if (pConfig->bDebug){
    Rcpp::Rcerr<<"Calling traverse events\n";
  }
  traverseEvents(true,xOverNode,newCoalEvent);
  if (pConfig->bDebug){
    Rcpp::Rcerr<<"Traverse events completed successfully\n";
  }
  dSplitPoint = newCoalEvent->getTime();
  bool bNewOrigin = false;
  EdgePtr selectedCoalEdge;
  if (dSplitPoint<this->grandMRCA->getHeight()){
    if(pConfig->bDebug) Rcpp::Rcerr<<"Coalescing at height "<<dSplitPoint<<endl;
    try{
      selectedCoalEdge = getRandomEdgeToCoalesce(
        coalescingEdge,dSplitPoint);
    }catch(unsigned int numEdges){
      printDataStructures();
      throw
        "Didn't find an edge to coalesce with, exiting\n";
    }
  }else{
    bNewOrigin = true;
  }
  NodePtr coalNode;
  // save the previous grand MRCA
  NodePtr oldGrandMRCA = grandMRCA;
  // discard the grandMRCA extending edge above the grandMRCA
  grandMRCA = originExtension->getBottomNodeRef();
  grandMRCA->clearTopEdges();
  if (bNewOrigin){
    coalNode = NodePtr(new Node(Node::COAL,
                                grandMRCA->getPopulation(),dSplitPoint));
    if(pConfig->bDebug) Rcpp::Rcerr<<"Creating new grand ancestor at "<<
      coalNode->getHeight()<<endl;
#ifdef DIAG
    if (grandMRCA->getPopulation()!=coalNode->getPopulation()){
      throw "Old grandMRCA and new grandMRCA don't match in pop\n";
    }
#endif
    EdgePtr shortEdge = EdgePtr(new Edge(coalNode,grandMRCA));
    shortEdge->iGraphIteration=iGraphIteration;
    addEdge(shortEdge);
    grandMRCA->addNewEdge(Node::TOP_EDGE,shortEdge);
    grandMRCA = coalNode;
    grandMRCA->addNewEdge(Node::BOTTOM_EDGE,shortEdge);
  }else{
    coalNode = NodePtr(new Node(
      Node::COAL,coalescingEdge->getBottomNodeRef()->getPopulation(),
      dSplitPoint));
    insertNodeInEdge(coalNode,selectedCoalEdge);
  }
  // set the variable proposed grandMRCA so that marking current tree can reach at least this
  // height
  localMRCA = coalNode;
  coalNode->setEvent(newCoalEvent);
  EdgePtr coalEdgeCopy = coalescingEdge;
  coalescingEdge = EdgePtr(new Edge(coalNode,coalescingEdge->getBottomNodeRef()));
  coalescingEdge->iGraphIteration = coalEdgeCopy->iGraphIteration;
  coalescingEdge->getBottomNodeRef()->replaceOldWithNewEdge(Node::TOP_EDGE,
                                   coalEdgeCopy,coalescingEdge);
  coalescingEdge->getTopNodeRef()->addNewEdge(Node::BOTTOM_EDGE,
                                coalescingEdge);
  addEdge(coalescingEdge);
  coalescingEdge->iGraphIteration=iGraphIteration;
  EdgePtr edgeBelowXover = xOverNode->getBottomEdgeByIndex(0);
  markEdgesBelow(edgeBelowXover);
}

void GraphBuilder::markCurrentTree(){
  iTotalTreeEdges = 0;
  unsigned int iTotalSamples=pConfig->iSampleSize;
  for(unsigned int i=0;i<iTotalSamples;++i){
    NodePtr & node = pSampleNodeArray[i];
    pTreeEdgesToCoalesceArray[i]=node->getTopEdgeByIndex(0);
  }
  unsigned int iIndex = 0;
  int iIterations = 0;
  bool bFirstSample = true;
  while(iIndex<iTotalSamples){
    EdgePtr & curEdge = pTreeEdgesToCoalesceArray[iIndex];
    if (markEdgesAbove(bFirstSample,true,curEdge,iIndex)){
      ++iIndex;
    }else{
      // start over
      iIndex = 0;
    }
    bFirstSample = false;
    ++iIterations;
  }
#ifdef DIAG
  if (localMRCA->getBottomEdgeByIndex(0)->iGraphIteration!=
      localMRCA->getBottomEdgeByIndex(1)->iGraphIteration){
    printDataStructures();
    Rcpp::Rcerr<<"Proposed grandMRCA at :"<<localMRCA->getHeight()<<endl;
    throw "Proposed grandMRCA's edges' histories don't match";
  }
#endif
}



void GraphBuilder::traverseEvents(bool bBuildFromEventList,
                                  NodePtr & xOverNode, EventPtr & newCoalEvent){
  short int iRequestedPop = -1;
  short int iOriginPop = -1;
  double dXoverHeight = 0.;
  if (bBuildFromEventList){
    iRequestedPop = xOverNode->getPopulation();
    iOriginPop = grandMRCA->getPopulation();
    dXoverHeight = xOverNode->getHeight();
  }
  // copy stuff over from Configuration object
  int iTotalPops = pConfig->iTotalPops;
  // make a copy of the population information to the dynamic vector
  PopVector pPopList = pConfig->pPopList;
  // make a copy of the the migration matrix over to our dynamic 2-d vector
  //if (pConfig->bMigrationChangeEventDefined){
  //dMigrationMatrix.clear();
  dMigrationMatrix = pConfig->dMigrationMatrix;
  //}
  if (pConfig->bDebug){
    //Rcpp::Rcerr<<"Migration Matrix from copy\n";
    //for (int j=0;j<iTotalPops;++j){
    //  for (int k=0;k<iTotalPops;++k){
    //    Rcpp::Rcerr<<" "<<dMigrationMatrix[j][j];
    //  }
    //  Rcpp::Rcerr<<endl;
    // }
  }
  // set up pile of coalesced nodes for building the prior tree
  NodePtrSet * pCoalescedNodes = NULL;
  if (!bBuildFromEventList){
    // store coalesced nodes in a temp vector
    //pCoalescedNodes = new NodePtrList();
    pCoalescedNodes = new NodePtrSet();
    // set up the node list
    int iCounter = 0,iId=0;
    for (int i=0;i<iTotalPops;++i){
      int iSampleSize = pConfig->pPopList[i].getChrSampled();
      for (int j=0;j<iSampleSize;++j){
        NodePtr newNode = NodePtr(new SampleNode
                                    (i,iId++));
        this->pSampleNodeArray[iCounter++]=newNode;
        pCoalescedNodes->insert(newNode);
      }
    }
  }
  
  // remove any leading events that are marked for deletion.
  if (pConfig->bDebug){
    Rcpp::Rcerr<<"Removing leading events marked for deletion\n";
  }
  EventPtrList::iterator currentEventIt = pEventList->begin();
  EventPtrList::iterator lastEventIt = pEventList->end();
  EventPtr pNextNewEvent;
  if (currentEventIt!=lastEventIt){
    pNextNewEvent = *currentEventIt;
    while (pNextNewEvent->bMarkedForDelete){
      currentEventIt = pEventList->erase(currentEventIt);
      pNextNewEvent = *currentEventIt;
    }
  }
  
  unsigned int iChromCounter = bBuildFromEventList?0:pCoalescedNodes->size();
  bool bDoneBuild = false,bIsEvent = false,
    bUserEventAvailable = false,bAboveOrigin = false,
    bAboveXoverHeight = false,bProposeNewCoal = false;
  short int iCoalescingPop = -1;
  unsigned long int iCoalRate = 0;
  double dTime = 0.0,dMigration = 0.0,dLastTime = 0.0,dWaitTime=0.;
  if (pConfig->bDebug){
    Rcpp::Rcerr<<"Beginning events traversal\n";
  }
  while (!bDoneBuild){
    int iEventType = -1;
    if (currentEventIt!=lastEventIt){
      pNextNewEvent = *currentEventIt;
      bUserEventAvailable = true;
      if (bBuildFromEventList){
        //Rcpp::Rcerr<<"Xover height "<<dXoverHeight<<endl;
        if (dXoverHeight>dLastTime &&
            dXoverHeight<=pNextNewEvent->getTime()){
          EventPtr eventWrapper = EventPtr(new XoverEvent(
            Event::XOVER,dXoverHeight,iRequestedPop));
          xOverNode->setEvent(eventWrapper);
          pEventList->insert(currentEventIt,eventWrapper);
          pNextNewEvent = eventWrapper;
          currentEventIt--;
        }
      }
    }else{
      bUserEventAvailable = false;
    }
    if (bBuildFromEventList){
      if (dTime==dXoverHeight){
        bAboveXoverHeight = true;
        bProposeNewCoal = true;
      }
      if (dTime==grandMRCA->getHeight()){
        bAboveOrigin = true;
      }
    }
    if (!bBuildFromEventList || bAboveXoverHeight){
      bIsEvent = false;
      double dProposedWaitTime = 0.;
      // BEGIN COALESCENT SECTION
      if(bProposeNewCoal ||!bBuildFromEventList){
        // attempt a coalescent event
        for(int pop=0; pop<iTotalPops ; ++pop) {      //coalescent
          int iChrSampled=(pPopList[pop].getPopSize()>0.)?pPopList[pop].getChrSampled():0;
          if (bProposeNewCoal){
            iCoalRate = iRequestedPop==pop?(iChrSampled-1)*2:0;
          }else{
            iCoalRate =  iChrSampled*(iChrSampled-1);
          }
          if( iCoalRate > 0 ) {
            if (pPopList[pop].getGrowthAlpha() == 0. ){
              dProposedWaitTime = pRandNumGenerator->expRV(
                iCoalRate/pPopList[pop].getPopSize());
            } else {  // support pop growth
              double dRandomVar=pRandNumGenerator->unifRV();
              double dCoalRateAlpha =
                1. - pPopList[pop].getGrowthAlpha()*
                pPopList[pop].getPopSize()*
                exp(-pPopList[pop].getGrowthAlpha()*
                (dTime - pPopList[pop].getLastTime() ) )*
                log(dRandomVar) / iCoalRate;
              if( dCoalRateAlpha > 0.0 ){
                // if  <= 0,  no coalescent within interval
                dProposedWaitTime = log(dCoalRateAlpha) /
                  pPopList[pop].getGrowthAlpha();
              }
            }
            if (!bIsEvent||dProposedWaitTime<dWaitTime) {
              iCoalescingPop = pop;
              dWaitTime = dProposedWaitTime;
            }
            iEventType = Event::NEW_COAL;
            bIsEvent = true;
            if (this->bEndGeneConversion){
              dWaitTime = 0.0;
            }
          }
        }
      }
      // END COALESCENT SECTION
      // BEGIN MIGRATION SECTION
      dMigration = 0.0;
      if (bBuildFromEventList){
        if (bAboveOrigin){
          dMigration = dMigrationMatrix[iRequestedPop][iRequestedPop]+
            dMigrationMatrix[iOriginPop][iOriginPop];
        }else{
          dMigration = dMigrationMatrix[iRequestedPop][iRequestedPop];
        }
      }else{
        for(int i=0; i<iTotalPops; ++i){
          dMigration +=
            pPopList[i].getChrSampled()*dMigrationMatrix[i][i];
        }
      }
      // Attempt a migration event
      if(dMigration > 0.0 ) {          // migration
        dProposedWaitTime = pRandNumGenerator->expRV(dMigration);
        if(!bIsEvent || dProposedWaitTime < dWaitTime){
          dWaitTime = dProposedWaitTime;
          iEventType = Event::NEW_MIGRATION;
          bIsEvent = true;
        }
      }else{
#ifdef DIAG
        if (dMigration<0.0)
          throw "Negative migration. Program exiting";
#endif
        if(iTotalPops>1 && !bUserEventAvailable){
          int iPops = 0;
          for(int j=0; j<iTotalPops; ++j) {
            if( pPopList[j].getChrSampled() > 0 ) ++iPops;
          }
          if( iPops > 1 ) {
            throw "Infinite coalescent time. No migration";
          }
        }
      }
      // END MIGRATION SECTION
    }
    // Now try one of the user invoked events
    if (bUserEventAvailable &&
        (!bIsEvent || dTime+dWaitTime > pNextNewEvent->getTime())
    ){
      dTime = pNextNewEvent->getTime();
      short int eventType =  pNextNewEvent->getType();
      //if (pConfig->bDebug) Rcpp::Rcerr<<"At time "<<dTime<<" event type is  "<<eventType<<endl;
      if (eventType==Event::PAST_COAL){
        CoalEvent * pExistingCoalEvent =
          static_cast<CoalEvent *>(pNextNewEvent.get());
        iCoalescingPop = pExistingCoalEvent->getPopulation();
        pPopList[iCoalescingPop].changeChrSampled(-1);
      }else if (eventType==Event::PAST_MIGRATION){
        MigrationEvent * pExistingMigEvent =
          static_cast<MigrationEvent *>(pNextNewEvent.get());
        short int iSourcePop = pExistingMigEvent->getPopMigratedFrom();
        short int iDestPop = pExistingMigEvent->getPopMigratedTo();
        
        // GKC 2015-07-05 disabled the code in condition below. 
        // Should not be necessary if POPSPLIT happens epsilon before
        // this time.
        iTotalPops = pPopList.size();
        if (false && iDestPop==iTotalPops){
          // we need to split the populations
          // GKC: 2015-07-01 Changed iTotalPops to iTotalPops-1
          vector <double> newRow(iTotalPops);
          for (int j=0;j<iTotalPops;++j){
            newRow[j]=pConfig->dGlobalMigration;
          }
          dMigrationMatrix.push_back(newRow);
          Population newPop;
          pPopList.push_back(newPop);
          iTotalPops = pPopList.size();
          for (int j=0;j<iTotalPops;++j){
            dMigrationMatrix[j].push_back(pConfig->dGlobalMigration);
          }
          if (pConfig->bDebug) Rcpp::Rcerr<<"dMigrationMatrix has dim "<<dMigrationMatrix.size()<<endl;
        }
        
        
        if ((iSourcePop>=iTotalPops||
            iDestPop>=iTotalPops)) {
          Rcpp::Rcerr<<"Invalid past migration event at time "
              <<dTime<<",source,dest pop "<<iSourcePop<<","<<
          iDestPop<<endl;
          throw "Num of pops is too small for migration!";
        }
        pPopList[iSourcePop].changeChrSampled(-1);
        pPopList[iDestPop].changeChrSampled(1);
        
      }else if (eventType==Event::XOVER){ //recombination
        XoverEvent * pXoverEvent =
          static_cast<XoverEvent *>(pNextNewEvent.get());
        pPopList[pXoverEvent->getPopulation()].changeChrSampled(1);
      }else if (eventType==Event::MIGRATION_MATRIX_RATE){
        MigrationRateMatrixEvent* migRateMatrixEvent = static_cast
        <MigrationRateMatrixEvent*>(pNextNewEvent.get());
        dMigrationMatrix = migRateMatrixEvent->getMigrationMatrix();
        if (dMigrationMatrix.size()!=pPopList.size()){
          Rcpp::Rcerr<<"Error in specifying new migration matrix event.  The dimension of the matrix "<<dMigrationMatrix.size()<<" must equal the number of populations: "<<pPopList.size()<<endl;
          throw "Invalid migration matrix";
        }
      }else if (eventType==Event::MIGRATION_RATE){
        MigrationRateEvent* migRateEvent = static_cast
        <MigrationRateEvent*>(pNextNewEvent.get());
        short int  iSourcePop = migRateEvent->getSourcePop();
        short int iDestPop = migRateEvent->getDestPop();
        double dMigRate = migRateEvent->getRate();
        if (pConfig->bDebug){
          Rcpp::Rcerr<<"Mig matrix has dim "<<dMigrationMatrix.size()<<" and source,Dest is "<<iSourcePop<<","<<iDestPop<<endl;
        }
        dMigrationMatrix[iSourcePop][iSourcePop]+=
          dMigRate-dMigrationMatrix[iSourcePop][iDestPop];
        dMigrationMatrix[iSourcePop][iDestPop] = dMigRate;
      }else if (eventType==Event::GLOBAL_MIGRATIONRATE){
        GenericEvent * pGeneriEvent =
          static_cast<GenericEvent *>(pNextNewEvent.get());
        for(short unsigned int i=0;i<iTotalPops;++i){
          for(short unsigned int j=0;j<iTotalPops;++j){
            dMigrationMatrix[i][j] =
              pGeneriEvent->getParamValue()/(iTotalPops-1);
          }
        }
        for(short unsigned int i=0;i<dMigrationMatrix.size();++i){
          dMigrationMatrix[i][i] =
            pGeneriEvent->getParamValue();
        }
      }else if (eventType==Event::GROWTH){
        PopSizeChangeEvent * growthEvent = static_cast
        <PopSizeChangeEvent *>(pNextNewEvent.get());
        short int pop = growthEvent->getPopulationIndex();
        pPopList[pop].setPopSize(pPopList[pop].getPopSize()*
          exp( - pPopList[pop].getGrowthAlpha()*(
              dTime-pPopList[pop].getLastTime()) ))  ;
        pPopList[pop].setGrowthAlpha(growthEvent->getPopChangeParam());
        pPopList[pop].setLastTime(dTime);
      }else if (eventType==Event::POPSIZE){
        PopSizeChangeEvent * growthEvent = static_cast
        <PopSizeChangeEvent *>(pNextNewEvent.get());
        short int pop = growthEvent->getPopulationIndex();
        pPopList[pop].setPopSize(growthEvent->getPopChangeParam());
        pPopList[pop].setGrowthAlpha(0);
      }else if (eventType==Event::POPJOIN){
        // merge pop i into pop j  (join)
        PopJoinEvent * joinEvent = static_cast
        <PopJoinEvent *>(pNextNewEvent.get());
        short int iSourcePop = joinEvent->getSourcePop();
        short int iDestPop = joinEvent->getDestPop();
        if (pConfig->bDebug){
          Rcpp::Rcerr<<"POPJOIN event encountered at time "<<dTime<<" from "<< iSourcePop<<" to "<<iDestPop <<"\n";
        }
        if (iSourcePop>=iTotalPops){
          Rcpp::Rcerr <<"Source pop and total pops are "<<iSourcePop<<","<<iTotalPops<<" in POP JOIN event at history "<<iGraphIteration<<". It is recommended that you increase the migration rates and/or number of sampled chromosomes.\n";
          
          throw "Invalid data structure";
        }
        
        if (iDestPop>=iTotalPops){
          Rcpp::Rcerr <<"Dest pop and total pops are "<<iDestPop<<","<<iTotalPops<<" in POP JOIN event\n";
          throw "Invalid data structure";
        }
        if (iDestPop==iTotalPops){
          // we need to split the populations
          Population newPop;
          pPopList.push_back(newPop);
          iTotalPops = pPopList.size();
          vector <double> newRow(iTotalPops);
          for (int j=0;j<iTotalPops;++j)
            newRow[j]=0.0;
          dMigrationMatrix.push_back(newRow);
          for (int j=0;j<iTotalPops;++j)
            dMigrationMatrix[j].push_back(0.0);
          if (pConfig->bDebug){
            Rcpp::Rcerr<<"In pop join, dMigrationMatrix dim is "<<dMigrationMatrix.size()<<endl;
          }
        }
        
        if (!bBuildFromEventList){
          NodePtrList nodesToMigrate;
          for (NodePtrSet::iterator it=pCoalescedNodes->begin();
               it!=pCoalescedNodes->end();++it){
            NodePtr node = *it;
            if (node->getPopulation()==iSourcePop){
              nodesToMigrate.push_back(node);
            }
          }
          for (NodePtrList::iterator it=nodesToMigrate.begin();
               it!=nodesToMigrate.end();++it){
            NodePtr node = *it;
            NodePtr childNode = node;
            NodePtr parentNode =
              NodePtr(new Node(Node::MIGRATION,
                               iDestPop,dTime));
            EdgePtr newEdge = EdgePtr(
              new Edge(parentNode,childNode));
            this->addEdge(newEdge);
            parentNode->addNewEdge(Node::BOTTOM_EDGE,newEdge);
            childNode->addNewEdge(Node::TOP_EDGE,newEdge);
            pCoalescedNodes->erase(childNode);
            pCoalescedNodes->insert(parentNode);
            // store this new event
            EventPtr eventWrapper = EventPtr(new MigrationEvent(
              Event::PAST_MIGRATION,dTime,
              iDestPop,iSourcePop));
            parentNode->setEvent(eventWrapper);
            pEventList->insert(currentEventIt,eventWrapper);
          }
        }else{
          if (iRequestedPop==iSourcePop &&
              coalescingEdge->getBottomNodeRef()->getHeight() < dTime &&
              dTime < coalescingEdge->getTopNodeRef()->getHeight()){
            NodePtr migrNode = NodePtr(new Node
                                         (Node::MIGRATION,iDestPop,dTime));
            insertNodeInRunningEdge(migrNode,coalescingEdge);
            iRequestedPop = iDestPop;
          }
          if (iOriginPop==iSourcePop &&
              originExtension->getBottomNodeRef()->getHeight() < dTime &&
              dTime < originExtension->getTopNodeRef()->getHeight()){
            NodePtr migrNode = NodePtr(new Node
                                         (Node::MIGRATION,iDestPop,dTime));
            insertNodeInRunningEdge(migrNode,originExtension);
            iOriginPop = iDestPop;
            // store this new event
            EventPtr eventWrapper = EventPtr(new MigrationEvent(
              Event::PAST_MIGRATION,dTime,
              iDestPop,iSourcePop));
            migrNode->setEvent(eventWrapper);
            pEventList->insert(currentEventIt,eventWrapper);
          }
        }
        pPopList[iDestPop].changeChrSampled(pPopList[iSourcePop]
                          .getChrSampled()) ;
        pPopList[iSourcePop].setChrSampled(0);
        // 2013-06-23 GKC begin codefix to revise migration matrix
        for (int j=0;j<iTotalPops;++j){
          if (j==iSourcePop){
            for (int k=0;k<iTotalPops;++k){
              // this population is dead and doesn't accept migrants
              dMigrationMatrix[j][k] = 0.;
            }
          }else{
            // look for the emigration rate of the dead pop
            for (int k=0;k<iTotalPops;++k){
              if (k==iSourcePop){
                dMigrationMatrix[j][j] -= dMigrationMatrix[j][k];
                dMigrationMatrix[j][k] = 0;;
                break;
              }
            }
          }
        }
        if (pConfig->bDebug){
          Rcpp::Rcerr<<"POPJOIN event completed\n";
        }
        // 2013-06-23 GKC end codefix to revise migration matrix
      }else if (eventType==Event::POPSPLIT){
        PopSizeChangeEvent *splitEvent = static_cast
        <PopSizeChangeEvent *>(pNextNewEvent.get());
        short int iSourcePop = splitEvent->getPopulationIndex();
        short int iDestPop = iTotalPops;
        if (pConfig->bDebug){
          Rcpp::Rcerr<<"POPSPLIT event encountered at time "<<dTime<<
            " from "<<iSourcePop<<" to "<<iDestPop<<"\n";
          //    if (iGraphIteration==887) throw "test exit";
        }
        double dProportion = splitEvent->getPopChangeParam();
        // GKC: 2015-07-05 Add epsilon to migration node
        // to make sure POPSPLIT happens first
        double dMigrationTime = dTime;
        // GKC: 2015-07-01 Changed iTotalPops to iTotalPops-1
        // begin expand migration matrix
        Population newPop;
        iTotalPops = pPopList.size();
        vector <double> newRow(iTotalPops);
        for (int j=0;j<(iTotalPops);++j)
          newRow[j]=pConfig->dGlobalMigration;
        dMigrationMatrix.push_back(newRow);
        pPopList.push_back(newPop);
        iTotalPops = pPopList.size();
        for (int j=0;j<iTotalPops;++j){
          dMigrationMatrix[j].push_back(pConfig->dGlobalMigration);
        }
        if(pConfig->bDebug) Rcpp::Rcerr<<"POPSPLIT migration matrix expanded "
                                <<"to "<<dMigrationMatrix.size()<<" rows.\n";
        // end expand migration matrix
        if (!bBuildFromEventList){
          NodePtrList nodesToMigrate;
          for (NodePtrSet::iterator it=pCoalescedNodes->begin();
               it!=pCoalescedNodes->end();++it){
            NodePtr node = *it;
            if (node->getPopulation()==iSourcePop &&
                pRandNumGenerator->unifRV()<dProportion){
              nodesToMigrate.push_back(node);
            }
          }
          //GKC 2015-07-05 make sure the migration event happens
          // after the pop split
          if(nodesToMigrate.size()>0) ++currentEventIt;
          for (NodePtrList::iterator it=nodesToMigrate.begin();
               it!=nodesToMigrate.end();++it){
            NodePtr node = *it;
            NodePtr childNode = node;
            NodePtr parentNode =
              NodePtr(new Node(Node::MIGRATION,
                               iDestPop,dMigrationTime));
            if(pConfig->bDebug) Rcpp::Rcerr<<
              "POPSPLIT in build tree, migration node at time "
              <<dMigrationTime<<
              endl;
            EdgePtr newEdge =
              EdgePtr(new Edge(parentNode,childNode));
            this->addEdge(newEdge);
            parentNode->addNewEdge(Node::BOTTOM_EDGE,newEdge);
            childNode->addNewEdge(Node::TOP_EDGE,newEdge);
            pCoalescedNodes->erase(childNode);
            pCoalescedNodes->insert(parentNode);
            pPopList[iDestPop].changeChrSampled(1) ;
            pPopList[iSourcePop].changeChrSampled(-1);
            // store this new event
            EventPtr eventWrapper = EventPtr(new MigrationEvent(
              Event::PAST_MIGRATION,dMigrationTime,
              iDestPop,iSourcePop));
            parentNode->setEvent(eventWrapper);
            pEventList->insert(currentEventIt,eventWrapper);
          }
          //GKC 2015-07-05 make sure the migration event happens
          // after the pop split
          if(nodesToMigrate.size()>0) --currentEventIt;
          if (pConfig->bDebug) {
            Rcpp::Rcerr<<"In pop split mig dim is size "<<dMigrationMatrix.size()<<endl;
          }
        }else{
          iDestPop=iTotalPops-1;
          if (pConfig->bDebug){
            Rcpp::Rcerr<<"In POPSPLIT: destpop,sourcepop,size: "<<iDestPop<<","<<iSourcePop<<","<<pPopList.size()<<endl;
          }
          double random = pRandNumGenerator->unifRV();
          bool internal_edge_condition =
            iRequestedPop==iSourcePop &&
            coalescingEdge->getBottomNodeRef()->getHeight() < 
              dMigrationTime && dMigrationTime < 
                coalescingEdge->getTopNodeRef()->getHeight();
          bool external_edge_condition = 
            iOriginPop&&
            originExtension->getBottomNodeRef()->getHeight() < 
              dMigrationTime && dMigrationTime < 
                originExtension->getTopNodeRef()->getHeight();
          
          if(random < dProportion && (internal_edge_condition ||
             external_edge_condition)){
            //if (iRequestedPop==iSourcePop&&
            //coalescingEdge->getBottomNodeRef()->getHeight() < 
            //dMigrationTime && dMigrationTime < coalescingEdge->getTopNodeRef()->getHeight()&&
            //    random<dProportion){
            NodePtr migrNode = NodePtr(new Node
                                         (Node::MIGRATION,iDestPop,dMigrationTime));
            if(internal_edge_condition){
              insertNodeInRunningEdge(migrNode,coalescingEdge);
              iRequestedPop = iDestPop;
            }
            if(external_edge_condition){
              insertNodeInRunningEdge(migrNode,originExtension);
              iOriginPop = iDestPop;
            }
            pPopList[iDestPop].changeChrSampled(1) ;
            pPopList[iSourcePop].changeChrSampled(-1);
            // store this new event
            EventPtr eventWrapper = EventPtr(new MigrationEvent(
              Event::PAST_MIGRATION,dMigrationTime,
              iDestPop,iSourcePop));
            migrNode->setEvent(eventWrapper);
            //GKC 2015-07-05 make sure the migration event happens
            // after the pop split
            ++currentEventIt;
            pEventList->insert(currentEventIt,eventWrapper);
            --currentEventIt;
          }
          //}
          //if (iOriginPop&&
          //originExtension->getBottomNodeRef()->getHeight() < 
          //dMigrationTime && dMigrationTime < 
          //originExtension->getTopNodeRef()->getHeight()&&
          //   random<dProportion){
          //NodePtr migrNode = NodePtr(new Node
          //(Node::MIGRATION,iDestPop,dMigrationTime));
          //insertNodeInRunningEdge(migrNode,originExtension);
          //iOriginPop = iDestPop;
          //pPopList[iDestPop].changeChrSampled(1) ;
          //pPopList[iSourcePop].changeChrSampled(-1);
          //// store this new event
          //EventPtr eventWrapper = EventPtr(new MigrationEvent(
          //Event::PAST_MIGRATION,dMigrationTime,
          //iDestPop,iSourcePop));
          //migrNode->setEvent(eventWrapper);
          ////GKC 2015-07-05 make sure the migration event happens
          //// after the pop split
          //++currentEventIt;
          //pEventList->insert(currentEventIt,eventWrapper);
          //}
        }
      }else if (eventType==Event::GLOBAL_POPGROWTH){
        GenericEvent * pGeneriEvent =
          static_cast<GenericEvent *>(pNextNewEvent.get());
        for(short unsigned int iPop=0;iPop<pPopList.size();++iPop){
          Population & pop = pPopList[iPop];
          pop.setPopSize(pop.getPopSize()*exp( - pop.getGrowthAlpha()*(dTime - pop.getLastTime())));
          pop.setGrowthAlpha(pGeneriEvent->getParamValue());
          pop.setLastTime(dTime);
        }
      }else if (eventType==Event::GLOBAL_POPSIZE){
        GenericEvent * pGeneriEvent =
          static_cast<GenericEvent *>(pNextNewEvent.get());
        for(short unsigned int iPop=0;iPop<pPopList.size();++iPop){
          Population & pop = pPopList[iPop];
          pop.setPopSize(pGeneriEvent->getParamValue());
          pop.setGrowthAlpha(0);
        }
      }else{
        Rcpp::Rcerr<<"Found an event of type "<<pNextNewEvent->getType()
            <<endl;
        throw "Event is not implemented yet!";
      }
      ++currentEventIt;
      //if (pConfig->bDebug) Rcpp::Rcerr<<"looking for events to mark for deletion\n";
      bool found = false;
      while(!found && currentEventIt!=lastEventIt){
        if (*currentEventIt==NULL) Rcpp::Rcerr<<"null\n";
        if ((*currentEventIt)->bMarkedForDelete) {
          if(pConfig->bDebug)Rcpp::Rcerr<<"Deleting event at "<<(*currentEventIt)->getTime()<<endl;
          currentEventIt=pEventList->erase(currentEventIt);
        }
        else found = true;
      }
      //if (pConfig->bDebug) Rcpp::Rcerr<<"user event processing complete\n";
    }else{
      dTime += dWaitTime;
      //if(pConfig->bDebug) Rcpp::Rcerr<<"Handling migration or coalescence at time "<<dTime<<endl;
      if(iEventType==Event::NEW_MIGRATION){
        double dRandMigr =  dMigration*pRandNumGenerator->unifRV();
        vector<int> migrantPops;
        if (!bBuildFromEventList){
          for (NodePtrSet::iterator it = pCoalescedNodes->begin();
               it!=pCoalescedNodes->end();++it){
            NodePtr node = *it;
            migrantPops.push_back(node->getPopulation());
          }
        }else{
          if (bAboveOrigin){
            migrantPops.push_back(int(iOriginPop));
          }
          // In any case a candidate migrant is the new coal line
          migrantPops.push_back(int(iRequestedPop));
        }
        int iTotalChrom = migrantPops.size();
        int migrant=-1,i=0;
        double dSum=0.;
        while (dRandMigr>=dSum && i<iTotalChrom){
          if (pConfig->bDebug){
            //Rcpp::Rcerr<<"In existing mig node, dMigrationMatrix size "<<dMigrationMatrix.size()<<" and migrantPops[i]: "<<migrantPops[i]<<endl;
          }
          dSum += dMigrationMatrix[migrantPops[i]][migrantPops[i]];
          migrant = i ;
          ++i;
        }
        int source_pop = -1,dest_pop = -1;
        NodePtr migrNode;
        if (!bBuildFromEventList){
          dRandMigr = pRandNumGenerator->unifRV()*
            dMigrationMatrix[migrantPops[migrant]][migrantPops[migrant]];
          dSum = 0.0;
          int i = 0;
          
          while (dRandMigr>=dSum && i<iTotalPops){
            if( i != migrantPops[migrant]){
              dSum += dMigrationMatrix[migrantPops[migrant]][i];
              dest_pop = i;
            }
            ++i;
          }
          NodePtrSet::iterator it = pCoalescedNodes->begin();
          for (int i=0;i<migrant;++i) ++it;
          NodePtr childNode = *it;
          source_pop = childNode->getPopulation();
          migrNode = NodePtr(new Node
                               (Node::MIGRATION,dest_pop,dTime));
          EdgePtr newEdge = EdgePtr(new Edge(migrNode,childNode));
          this->addEdge(newEdge);
          migrNode->addNewEdge(Node::BOTTOM_EDGE,newEdge);
          childNode->addNewEdge(Node::TOP_EDGE,newEdge);
          pCoalescedNodes->erase(childNode);
          pCoalescedNodes->insert(migrNode);
        }else{
          if (bAboveOrigin){
            source_pop = migrantPops[migrant];
          }else{
            source_pop=iRequestedPop;
          }
          
          dRandMigr = pRandNumGenerator->unifRV()*
            dMigrationMatrix[source_pop][source_pop];
          dSum = 0.0;
          i = 0;
          if (pConfig->bDebug){
            Rcpp::Rcerr<<"2 In existing mig node, dMigrationMatrix size "<<dMigrationMatrix.size()<<" and sourcepop, itotalpops: "<<source_pop<<","<<iTotalPops<<endl;
          }
          
          while (dRandMigr>=dSum && i<iTotalPops){
            if( i != source_pop){
              dSum += dMigrationMatrix[source_pop][i];
              dest_pop = i;
            }
            ++i;
          }
          migrNode = NodePtr(new Node
                               (Node::MIGRATION,dest_pop,dTime));
          if (bAboveOrigin){
            // There are two lines above the grandMRCA and either one
            // can have a migration event inserted.
            switch(migrant){
            case 0:
              insertNodeInRunningEdge(migrNode,originExtension);
              iOriginPop = dest_pop;
              break;
            case 1:
              insertNodeInRunningEdge(migrNode,coalescingEdge);
              iRequestedPop = dest_pop;
              break;
            default:
              throw "The index of the migration source is out of bounds.";
            }
          }else{
            insertNodeInRunningEdge(migrNode,coalescingEdge);
            iRequestedPop = dest_pop;
          }
        }
        pPopList[source_pop].changeChrSampled(-1);
        pPopList[dest_pop].changeChrSampled(1);
        EventPtr eventWrapper = EventPtr(new MigrationEvent(
          Event::PAST_MIGRATION,dTime,
          dest_pop,source_pop));
        migrNode->setEvent(eventWrapper);
        pEventList->insert(currentEventIt,eventWrapper);
      }else if (iEventType==Event::NEW_COAL){ // coalescent event
        
        if (!bBuildFromEventList){
          EventPtr eventWrapper;
          // pick the two, c1, c2
          int iChr1,iChr2;
          NodePtr node0,node1;
          bool found = false;
          NodePtr chr1Node,chr2Node;
          
          while(!found){
            iChr1 = static_cast<int>(pRandNumGenerator->unifRV() *
              pCoalescedNodes->size());
            NodePtrSet::iterator it = pCoalescedNodes->begin();
            for (int i=1;i<=iChr1;++i) ++it;
            chr1Node = *it;
            if (chr1Node->getPopulation()==iCoalescingPop){
              while(!found){
                iChr2 = static_cast<int>(pRandNumGenerator->unifRV() *
                  pCoalescedNodes->size());
                if (iChr1!=iChr2){
                  it = pCoalescedNodes->begin();
                  for (int i=1;i<=iChr2;++i) ++it;
                  chr2Node = *it;
                  if (chr2Node->getPopulation()==iCoalescingPop) found = true;
                }
              }
            }
          }
          node0 = chr1Node;
          node1 = chr2Node;
          NodePtr parentNode = NodePtr(new Node
                                         (Node::COAL,node0->getPopulation(),dTime));
          EdgePtr edge0 = EdgePtr(new Edge(parentNode,node0));
          node0->addNewEdge(Node::TOP_EDGE,edge0);
          parentNode->addNewEdge(Node::BOTTOM_EDGE,edge0);
          EdgePtr edge1 = EdgePtr(new Edge(parentNode,node1));
          node1->addNewEdge(Node::TOP_EDGE,edge1);
          parentNode->addNewEdge(Node::BOTTOM_EDGE,edge1);
          this->addEdge(edge0);
          this->addEdge(edge1);
          pCoalescedNodes->erase(node0);
          pCoalescedNodes->erase(node1);
          eventWrapper = EventPtr(new CoalEvent(
            Event::PAST_COAL,dTime,
            parentNode->getPopulation()));
          pCoalescedNodes->insert(parentNode);
          iChromCounter = pCoalescedNodes->size();
          if (iChromCounter==1){
            localMRCA = grandMRCA = parentNode;
            bDoneBuild = true;
          }
          pPopList[iCoalescingPop].changeChrSampled(-1);
          parentNode->setEvent(eventWrapper);
          pEventList->insert(currentEventIt,eventWrapper);
        }else{
          // the case where the proposed time is beyond the last
          // event time
          // or the proposed time is in the middle and an eligible
          // coalescent event is found.
          newCoalEvent = EventPtr(new CoalEvent(
            Event::PAST_COAL,dTime,iRequestedPop));
          pEventList->insert(currentEventIt,newCoalEvent);
          bDoneBuild = true;
        }
      }else{
        if (pConfig->bDebug) Rcpp::Rcerr<< "No random event could be assigned.";
      }
      
    }
    dLastTime = dTime;
  }
  if (!bBuildFromEventList){
    delete pCoalescedNodes;
  }
}

void GraphBuilder::pruneARG(int iHistoryMax){
  bool found = false;
  EdgePtrList candidateEdges;
  EdgePtrList::iterator it = pEdgeListInARG->begin();
  // look for old edges that have an expired graph iteration value
  while (it!=pEdgeListInARG->end()){
    EdgePtr curEdge = *it;
    if (curEdge->iGraphIteration<=iHistoryMax)
      candidateEdges.push_back(curEdge);
    ++it;
  }
  candidateEdges.sort(byBottomNode());
  do{
    EdgePtr oldEdge;
    found = false;
    // sort the edges from top to bottom as the pruning order
    EdgePtrList::iterator it2 = candidateEdges.begin();
    while(!found && it2!=candidateEdges.end()){
      EdgePtr curEdge = *it2;
      if (!curEdge->bDeleted &&curEdge->getBottomNodeRef()->
          getType()==Node::XOVER){
        oldEdge = curEdge;
        found = true;
      }else{
        ++it2;
      }
    }
    if (found){
      NodePtr & xoverNode = oldEdge->getBottomNodeRef();
      EdgePtr edge0 = xoverNode->getTopEdgeByIndex(0);
      EdgePtr edge1 = xoverNode->getTopEdgeByIndex(1);
      EdgePtr edgeToPrune = edge0->iGraphIteration<
        edge1->iGraphIteration?edge0:edge1;
      
      pruneEdgesAbove(edgeToPrune);
      
      edge0 = xoverNode->getTopEdgeByIndex(0);
      edge1 = xoverNode->getTopEdgeByIndex(1);
      
      EdgePtr otherTopEdge = edge0->iGraphIteration>
        edge1->iGraphIteration?edge0:edge1;
      
      EdgePtr bottomEdge = xoverNode->getBottomEdgeByIndex(0);
      
      if (grandMRCA->getHeight()>xoverNode->getHeight()){
        mergeEdges(otherTopEdge,bottomEdge);
        xoverNode->removeEvent();
      }else{
        // we just removed a loop above the grandMRCA which
        // includes the edges around the xover point.
      }
    }
  }while(found);
}

void GraphBuilder::addMutations(double startPos,double endPos){
  bool bEndMutate = false;
  while(!bEndMutate){
    // find the next point on this interval
    startPos+=pRandNumGenerator->expRV(dLastTreeLength*
      pConfig->dTheta);
    if (startPos>=endPos){
      bEndMutate = true;
    }else{
      double dRandomSpot = pRandNumGenerator->unifRV() * dLastTreeLength;
      double dMutationTime=-1.;
      EdgePtr selectedEdge = getRandomEdgeOnTree(dMutationTime,dRandomSpot);
      //Rcpp::Rcerr<<"Mutation time is "<<dMutationTime<<endl;
      mutateBelowEdge(selectedEdge);
      NodePtrVector::iterator it;
      
      
#ifdef DIAG
      Rcpp::Rcout<<MUTATIONSITE<<FIELD_DELIMITER<<pMutationPtrVector->size()<<
        FIELD_DELIMITER<< setw(15) << setprecision(9)<<startPos<<
          FIELD_DELIMITER<< setw(15) << setprecision(9)<<dMutationTime<<
            FIELD_DELIMITER;
#endif
      
      unique_ptr<AlphaSimRReturn> temp(new AlphaSimRReturn());
      temp->length = startPos;
      unsigned int iSampleSize = pConfig->iSampleSize;
      for (unsigned int iSampleIndex=0;iSampleIndex<iSampleSize;++iSampleIndex){
        SampleNode * sample = static_cast<SampleNode*>(pSampleNodeArray[iSampleIndex].get());
        sites[iSampleIndex]=sample->bAffected;
#ifdef DIAG
        Rcpp::Rcout<<sample->bAffected;
#endif
        temp->haplotypes.push_back(sample->bAffected);
        sample->bAffected=false;
      }
#ifdef DIAG
      Rcpp::Rcout<<endl;
#endif
      mutations.push_back(*temp);
      double dFreq=0.;
      if (pConfig->bSNPAscertainment){
        // first compute the MAF
        int counts=0;
        for (unsigned int i=0;i<iSampleSize;++i){
          counts=counts+sites[i];
        }
        dFreq = 1.*counts/iSampleSize;
        if (pConfig->bFlipAlleles && dFreq>.5){
          for (unsigned int i=0;i<iSampleSize;++i){
            sites[i]=!sites[i];
          }
          dFreq = 1.-dFreq;
        }
        AlleleFreqBinPtr query(new AlleleFreqBin(dFreq,dFreq,0.));
        AlleleFreqBinPtrSet::iterator it = pConfig->pAlleleFreqBinPtrSet->find(query);
        if (it!=pConfig->pAlleleFreqBinPtrSet->end()){
          AlleleFreqBinPtr bin = *it;
          ++bin->iObservedCounts;
        }else throw "Did not find a frequency range for freq";
      }
      pMutationPtrVector->push_back(new Mutation(startPos, dFreq));
      
    }
  }
}

bool GraphBuilder::getNextPos(double & curPos,HotSpotBinPtrList::iterator & hotSpotIt){
  bool bBinCrossed = false;
  if (hotSpotIt==pConfig->pHotSpotBinPtrList->end()){
    curPos+=pRandNumGenerator->expRV(getRate());
  }else{
    HotSpotBinPtr hotspot = *hotSpotIt;
    double startPos = hotspot->dStart;
    double endPos = hotspot->dEnd;
    if (curPos>=startPos&&curPos<endPos){
      curPos+=pRandNumGenerator->expRV(getRate()*hotspot->dRate);
      if (curPos>endPos){
        bBinCrossed = true;
        curPos = endPos;
        ++hotSpotIt;
      }
    }else if (curPos<startPos){
      curPos+=pRandNumGenerator->expRV(getRate());
      if (curPos>startPos){
        bBinCrossed = true;
        curPos = startPos;
      }
    }else{
      Rcpp::Rcerr<<"startPos is "<<startPos<<" endPos is "
          <<endPos<<" and curPos is "<<curPos<<endl;
      throw "Shouldn't be here for variable recomb";
    }
  }
  return bBinCrossed;
}

bool GraphBuilder::checkPendingGeneConversions(double & curPos){
  bool found = false;
  GeneConversionPtrSet::iterator it = pGeneConversionPtrSet->begin();
  while(!found){
    if (it!=pGeneConversionPtrSet->end()){
      GeneConversionPtr gc= *it;
      if (gc->dEndPos<curPos){
        EdgePtr edge0 = gc->xOverNode->getTopEdgeByIndex(0);
        EdgePtr edge1 = gc->xOverNode->getTopEdgeByIndex(1);
        if (gc->xOverNode->bDeleted ||
            edge1->iGraphIteration!=iGraphIteration){
          if (pConfig->bDebug){
            Rcpp::Rcerr<<"Deleting gene conversion because:"<<endl;
            if (gc->xOverNode->bDeleted)
              Rcpp::Rcerr<<"xover node was deleted\n";
            else if (edge1->iGraphIteration!=iGraphIteration)
              Rcpp::Rcerr<<"new edge "<<edge1->getBottomNodeRef()->getHeight()<<" to "<<
                edge1->getTopNodeRef()->getHeight()<<" was not ancestral\n";
          }
          delete(*it);
          pGeneConversionPtrSet->erase(it++);
        }else{
          if (pConfig->bDebug){
            Rcpp::Rcerr<<"Closing a gene conversion: "<<
              "Backtracking position from "<<curPos<<" to "<<
                gc->dEndPos<<endl;
          }
          curPos = gc->dEndPos;
          this->gcOldEdge = edge0;
          this->gcNewEdge = edge1;
          delete(*it);
          pGeneConversionPtrSet->erase(it++);
          return true;
        }
      }else{
        found = true;
      }
    }else found = true;
  }
  return false;
}

void GraphBuilder::build(){
  double curPos = 0.0,lastPos = 0.0,dMaxPos = 1.0;
  unsigned int iLastCumulativePos = 0;
#ifdef DIAG
  Rcpp::Rcerr<<"Debugging: "<<pConfig->bDebug<<endl;
#endif
  
  HotSpotBinPtrList::iterator hotSpotIt;
  if (pConfig->bVariableRecomb){
    hotSpotIt=pConfig->pHotSpotBinPtrList->begin();
  }
  // gene conversion stuff
  GeneConversionPtr newGC;
  double dLogTractRatio = log((pConfig->iGeneConvTract-1.)/pConfig->iGeneConvTract);
  const int XOVER=0,GCSTART=1,GCEND=2;
  int xovers=0,gcstarts=0,gcends=0;
  int lastRE;
  int iHistoryMax = 0;
  do{
    if (iGraphIteration==0){
      NodePtr dummy1;
      EventPtr dummy2;
      this->traverseEvents(false,dummy1,dummy2);
    }else{
      // at this point decide whether we invoke a plain x-over
      // or a new gene conversion event
      this->bBeginGeneConversion = false;
      lastRE=XOVER;
      if (this->bEndGeneConversion){
        lastRE=GCEND;
      }else{
        this->bBeginGeneConversion = pRandNumGenerator->unifRV()<
          (pConfig->dGeneConvRatio/(pConfig->dGeneConvRatio+1))?true:false;
        if (bBeginGeneConversion){
          lastRE=GCSTART;
          double dTractLen = (1.+log(pRandNumGenerator->unifRV())/
                              dLogTractRatio)/pConfig->dSeqLength;
          if (pConfig->bDebug) Rcpp::Rcerr<<"Proposing tract length: "
                                   <<dTractLen<<endl;
          newGC = GeneConversionPtr(new GeneConversion(
            curPos+dTractLen));
          pGeneConversionPtrSet->insert(newGC);
        }
      }
      invokeRecombination(newGC);
      switch(lastRE){
      case XOVER:
        ++xovers;
        break;
      case GCSTART:
        ++gcstarts;
        break;
      case GCEND:
        ++gcends;
        break;
      }
      // mark the graph edges as the local tree
      markCurrentTree();
      if (!bIncrementHistory){
        double dBoundary = curPos - dTrailingGap;
        if (dBoundary>0.){
          bIncrementHistory = true;
        }
      }else{
        ++iHistoryMax;
      }
#ifdef DIAG
      Rcpp::Rcerr<<"iHistoryMax: "<<iHistoryMax<<endl;
#endif
      if (iHistoryMax>=0){
        pruneARG(iHistoryMax);
      }
    }
    
    initializeCurrentTree();
#ifdef DIAG
    Rcpp::Rcerr<<"Tree:"<<iGraphIteration<<
      ",pos:"<<curPos<<
        ",len:"<<dLastTreeLength<<
          ",TMRCA:"<<localMRCA->getHeight()<<
            ",ARG:"<<
              ",len:"<<dArgLength<<
                ",TMRCA:"<<grandMRCA->getHeight()<<
                  endl;
#endif
    if (pConfig->bVariableRecomb){
      bool bBinCrossed;
      do{
        bBinCrossed = getNextPos(curPos,hotSpotIt);
      }while(bBinCrossed);
    }else{
      curPos+=pRandNumGenerator->expRV(getRate());
    }
    // check if we reached the end of the region
    if (curPos>dMaxPos) curPos=dMaxPos;
    if (pConfig->bNewickFormat){
      // Chen
      //uint iSegLength = (curPos-lastPos)*pConfig->dSeqLength;
      // Schiffels following two lines fixes
      uint iSegLength = curPos*pConfig->dSeqLength-iLastCumulativePos;
      iLastCumulativePos += iSegLength;
#ifdef DIAG
      Rcpp::Rcout<<NEWICKTREE<<"\t["<<iSegLength<<"]"<<
        getNewickTree(localMRCA->getHeight(),localMRCA)<<";"<<endl;
#endif
    }
    // check if there was an existing gene conversion event that needs
    // to be closed. backtrack if necessary.
    this->bEndGeneConversion  = checkPendingGeneConversions(curPos);
    if (pConfig->dTheta>0.0){
      addMutations(lastPos,curPos);
    }
    lastPos = curPos;
    ++iGraphIteration;
  }while(curPos<dMaxPos);
#ifdef DIAG
  Rcpp::Rcerr<<"Completed the chromosome at position "<<curPos<<endl;
  //if (pMutationPtrVector->size()>0){
  //        Rcpp::Rcout<<HAPLOEND<<endl;
  //}
  Rcpp::Rcerr<<"gcstarts:"<<gcstarts<<" gcends:"<<gcends<<" xovers:"<<xovers<<endl;
#endif
}

vector<AlphaSimRReturn> GraphBuilder::getMutations() {
  return mutations;
}

