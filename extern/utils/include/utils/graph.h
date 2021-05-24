/**
 * @file graph.h
 * @brief Basic graph data structure
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2017
 */

#include "asserter.h"
#include "pq.h"
#include <list>
#include <set>
#include <vector>
#include <limits>

#ifndef GRAPH_H
#define GRAPH_H

namespace dominiqs
{

template<typename EdgeData>
double get_weight(const EdgeData& data)
{
	return data;
}

/**
 * Undirected vertex-colored graph
 */
template <typename NodeData, typename EdgeData>
class Graph
{
public:
	struct Edge
	{
	public:
		Edge() = default;
		Edge(int s, int t, const EdgeData& d) : source(s), target(t), data(d) {}
		int source = 0;
		int target = 0;
		EdgeData data;
	};
	/**
	 * Add nodes/edges
	 */
	void addNode(NodeData n)
	{
		nodes.emplace_back(n);
	}
	void addEdge(int source, int target, EdgeData e)
	{
		// NOTE: we support multiple edges for the same pair of nodes
		edges.emplace_back(source, target, e);
	}
	void addUndirectedEdge(int source, int target, EdgeData e)
	{
		// NOTE: we support multiple edges for the same pair of nodes
		edges.emplace_back(source, target, e);
		edges.emplace_back(target, source, e);
	}
	/**
	 * Digest edges so that they are sorted by (source,target).
	 * start is initialized by this method.
	 */
	void digest()
	{
		// init start
		start.resize(nodes.size()+1);
		std::fill(start.begin(), start.end(), 0);
		// init adjacency lists
		std::vector<Edge> adjList(edges.size());
		// bucket sort
		std::vector<int> counts(nodes.size(), 0);
		for (const Edge& e: edges) counts[e.source]++;
		// compute starts
		int offset = 0;
		for (size_t v = 0; v < nodes.size(); v++)
		{
			start[v] = offset;
			offset += counts[v];
		}
		start[nodes.size()] = edges.size();
		// put elements into place (fill bucket backwards)
		for (const Edge& e: edges)
		{
			int node = e.source;
			adjList[start[node] + (--counts[node])] = e;
		}
		edges.swap(adjList);
	}
	/**
	 * Data
	 *
	 * The star of node i is made by the edges in position
	 * start[i]..start[i+1].
	 */
	std::vector<NodeData> nodes;
	std::vector<size_t> start;
	std::vector<Edge> edges;
	/**
	 * Dijkstra's shortest path algorithm
	 *
	 * Users need to implement the function get_weight(const EdgeData& d)
	 * for the algorithm to be able to extract the distance from each edge.
	 */
	void shortestPath(int source, std::vector<int>& prec, std::vector<double>& dist)
	{
		// initialize prec and dist
		size_t numNodes = nodes.size();
		prec.resize(numNodes);
		dist.resize(numNodes);
		std::fill(prec.begin(), prec.end(), -1);
		std::fill(dist.begin(), dist.end(), std::numeric_limits<double>::max());
		// add source to queue
		PriorityQueue<double> pq(numNodes);
		dist[source] = 0.0;
		pq.push(source, 0.0);
		// main loop
		while (!pq.isEmpty())
		{
			int u = pq.top();
			pq.pop();
			// traverse all edges in the star of u
			for (size_t itr = start[u]; itr < start[u+1]; itr++)
			{
				const Edge& edge = edges[itr];
				DOMINIQS_ASSERT(edge.source == u);
				int v = edge.target;
				double w_uv = get_weight(edge.data);
				if ( dist[u] + w_uv < dist[v] )
				{
					dist[v] = dist[u] + w_uv;
					prec[v] = u;
					if (pq.has(v))
					{
						// node already in the queue: check if we can decrease its distance
						pq.change(v, dist[v]);
					}
					else
					{
						// new node
						pq.push(v, dist[v]);
					}
				}
			}
		}
	}
	/**
	 * Bellmann-Ford's shortest path algorithm
	 *
	 * Users need to implement the function get_weight(const EdgeData& d)
	 * for the algorithm to be able to extract the distance from each edge.
	 */
	void shortestPathBellmannFord(int source, std::vector<int>& prec, std::vector<double>& dist, std::vector<size_t>& hops, size_t maxHop = -1)
	{
		// initialize prec and dist
		double infinity = std::numeric_limits<double>::max();
		size_t numNodes = nodes.size();
		prec.resize(numNodes);
		dist.resize(numNodes);
		hops.resize(numNodes);
		std::fill(prec.begin(), prec.end(), -1);
		std::fill(dist.begin(), dist.end(), infinity);
		std::fill(hops.begin(), hops.end(), -1);
		// set source to queue
		dist[source] = 0.0;
		hops[source] = 0;
		// hops loop
		for (size_t hop = 0; hop < maxHop; hop++)
		{
			// edges loop: we just need to loop over all the edges outgoing from a node currently reachable from the source
			for (size_t u = 0; u < numNodes; u++)
			{
				if (dist[u] == infinity)  continue;
				for (size_t itr = start[u]; itr < start[u+1]; itr++)
				{
					const Edge& edge = edges[itr];
					DOMINIQS_ASSERT(edge.source == u);
					int v = edge.target;
					double w_uv = get_weight(edge.data);
					if ( dist[u] + w_uv < dist[v] )
					{
						dist[v] = dist[u] + w_uv;
						prec[v] = u;
					}
				}
			}
		}
		// compute # of hops in shortest paths (basically a breadth-first)
		std::vector<int> queue(numNodes);
		int begin = 0;
		int end = 0;
		queue[end++] = source;
		while (begin < end)
		{
			int u = queue[begin++];
			if (dist[u] == infinity)  continue;
			for (size_t itr = start[u]; itr < start[u+1]; itr++)
			{
				const Edge& edge = edges[itr];
				int v = edge.target;
				if (prec[v] == u)
				{
					DOMINIQS_ASSERT(hops[v] == -1);
					hops[v] = hops[u] + 1;
					queue[end++] = v;
				}
			}
		}
	}
};

} // namespace dominiqs

#endif /* end of include guard: GRAPH_H */
