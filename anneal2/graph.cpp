#include "graph.h"

Node::~Node()
{
    EdgeList::iterator i = edges.begin();
    while(i != edges.end())
        delete *i++;
}

Edge::Edge(Node *a, Node *b) : nodeA(a), nodeB(b)
{
    iterA = nodeA->edges.insert(nodeA->edges.end(), const_cast<Edge *>(this));
    iterB = nodeB->edges.insert(nodeB->edges.end(), this);
}

Edge::~Edge()
{
    nodeA->edges.erase(iterA);
    nodeB->edges.erase(iterB);
}
