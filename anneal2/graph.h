#ifndef GRAPH_H
#define GRAPH_H

#include <list>

using std::list;


class Edge;

class Node
{
    friend class Edge;
    typedef list<Edge *> EdgeList;

    public:
        ~Node();

    private:
        Node(const Node &); // Disabled.
        Node &operator = (const Node &);  // Disabled.

    private:
        list<Edge *> edges;
};

class Edge
{
    public:
        Edge(Node *a, Node *b);
        ~Edge();

    private:
        Node *nodeA, *nodeB;
        Node::EdgeList::iterator iterA, iterB;
};

#endif
