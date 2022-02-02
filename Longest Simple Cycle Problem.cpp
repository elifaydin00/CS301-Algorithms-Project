#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>

using namespace std;

template <class obj> class StackNode;
template <class obj> class Stack;
enum class Color;
struct Vertex;
void LSC_BRUTEFORCE(int, int, vector<int>, const vector<vector<int>> &, vector<Stack<int>> &, Stack<int>);
void LSC(const int &, vector<Vertex> &, const vector<vector<int>> &);
void LSC_TRAV(const int &, const int &, vector<Vertex> &, const vector<vector<int>> &, int);
void PRINT_CYCLE(const int &, vector<Vertex> &);

template <class obj>
class StackNode
{
public:
	StackNode() : next(NULL) {}
	StackNode(const obj &_object, StackNode<obj> *_next) : object(_object), next(_next) {}
	obj object;
	StackNode<obj> *next;
};

template <class obj>
class Stack
{
public:
	Stack() : size(0), top(NULL) {}
	Stack(const Stack<obj> &rhs)
	{
		top = NULL;
		*this = rhs;
	}
	~Stack()
	{
		makeEmpty();
	}
	bool isEmpty() const
	{
		return top == NULL;
	}
	int getSize() const
	{
		return size;
	}
	void makeEmpty()
	{
		while (!isEmpty())
		{
			pop();
		}
		size = 0;
	}
	const obj &getTop() const
	{
		return top->object;
	}
	void push(const obj &_object)
	{
		top = new StackNode<obj>(_object, top);
		size++;
	}
	const obj &pop()
	{
		StackNode<obj> *prevTop = top;
		top = top->next;
		obj object = prevTop->object;
		delete prevTop;
		size--;
		return object;
	}
	const Stack<obj> &operator=(const Stack<obj> &rhs)
	{
		if (this != &rhs)
		{
			makeEmpty();
			if (rhs.isEmpty())
			{
				return *this;
			}
			Stack<obj> dummyStack;
			for (StackNode<obj> *rptr = rhs.top; rptr != NULL; rptr = rptr->next)
			{
				dummyStack.push(rptr->object);
			}
			while (!dummyStack.isEmpty())
			{
				obj object = dummyStack.pop();
				push(object);
			}
			size = rhs.size;
		}
		return *this;
	}

private:
	int size;
	StackNode<obj> *top;
};

enum class Color { White, Pink, Red };

struct Vertex
{
	Color color;
	int discover;
	int finish;
	int pred;
	int vertex_code;
	int adj_counter;
	Vertex() : color(Color::White), discover(-1), finish(-1), pred(-1), vertex_code(-1), adj_counter(-1) {}
	Vertex(int _vertex_code, int _adj_counter) : color(Color::White), discover(-1), finish(-1), pred(-1), vertex_code(_vertex_code), adj_counter(_adj_counter) {}
};

void LSC_BRUTEFORCE(int start, int last, vector<int> isVisited, const vector<vector<int>> &adjacencyMatrix, vector<Stack<int>> &paths, Stack<int> currentStack)
{
	int current = currentStack.getTop();
	if (current == start && isVisited[current] == 1)
	{
		return;
	}
	isVisited[current] = 1;
	for (int index = 0; index < (int)adjacencyMatrix[current].size(); index++)
	{
		if (adjacencyMatrix[current][index] == 1 && (isVisited[index] == 0 || index == start) && last != index)
		{
			currentStack.push(index);
			LSC_BRUTEFORCE(start, current, isVisited, adjacencyMatrix, paths, currentStack);
			if (currentStack.getTop() == start)
			{
				paths.push_back(currentStack);
			}
			currentStack.pop();
		}
	}
}

void LSC(const int &root, vector<Vertex> &vertices, const vector<vector<int>> &adj_matrix)
{
	for (int i = 0; i < (int)vertices.size(); i++)
	{
		if (vertices[i].vertex_code != root)
		{
			vertices[i].color = Color::White;
		}
		else
		{
			vertices[i].color = Color::Pink;
		}
		vertices[i].pred = -1;
	}
	int count = 0;
	for (int i = 0; i < (int)vertices.size(); i++)
	{
		if (adj_matrix[i][root] == 1 && vertices[i].color == Color::White)
		{
			vertices[i].pred = root;
			LSC_TRAV(i, root, vertices, adj_matrix, count);
		}
	}
}

void LSC_TRAV(const int &current, const int &root, vector<Vertex> &vertices, const vector<vector<int>> &adj_matrix, int count)
{
	vertices[current].color = Color::Pink;
	count++;
	vertices[current].discover = count;
	for (int i = 0; i < (int)vertices.size(); i++)
	{
		if (adj_matrix[i][current] == 1)
		{
			if (i == root)
			{
				vertices[i].finish = count + 1;
			}
			else if (vertices[i].adj_counter == 1)
			{
				vertices[i].finish = count;
			}
			else if (vertices[i].color == Color::White)
			{
				vertices[i].pred = current;
				LSC_TRAV(i, root, vertices, adj_matrix, count);
			}
		}
	}
	vertices[current].color = Color::Red;
	count++;
	vertices[current].finish = count;
}

void PRINT_CYCLE(const int &root, vector<Vertex> &vertices)
{
	int finish_max_index = 0;
	for (int i = 1; i < (int)vertices.size(); i++)
	{
		if (vertices[i].finish > vertices[finish_max_index].finish)
		{
			finish_max_index = i;
		}
	}
	int itr = finish_max_index;

	cout << root;
	while (vertices[itr].pred != -1)
	{
		cout << " " << itr;
		itr = vertices[itr].pred;
	}
	cout << endl;
}

float STD_DEV_CALCULATE(vector<float> &data)
{
	float total = 0.0;
	float mean;
	float stdDev = 0.0;
	for (int i = 0; i < (int)data.size(); i++)
	{
		total += data[i];
	}
	mean = total / (int)data.size();
	for (int i = 0; i < (int)data.size(); i++)
	{
		stdDev += pow(data[i] - mean, 2);
	}
	stdDev = sqrt(stdDev / (int)data.size());
	return stdDev;
}

float STD_ERR_CALCULATE(float stdDev, int size)
{
	return stdDev / sqrt(size);
}

void STATS_DISPLAY(vector<float> &rt)
{
	float total = 0.0;
	for (int i = 0; i < (int)rt.size(); i++)
	{
		total += rt[i];
	}

	float mean = total / (int)rt.size();
	float stdDev = STD_DEV_CALCULATE(rt);
	float stdErr = STD_ERR_CALCULATE(stdDev, (int)rt.size());

	const float tval90 = 1.645;
	float upperMean90 = mean + (tval90 * stdErr);
	float lowerMean90 = mean - (tval90 * stdErr);

	const float tval95 = 1.96;
	float upperMean95 = mean + (tval95 * stdErr);
	float lowerMean95 = mean - (tval95 * stdErr);

	cout << fixed << setprecision(6) << mean << " "
		<< fixed << setprecision(6) << stdDev << " "
		<< fixed << setprecision(6) << stdErr << " "
		<< fixed << setprecision(6) << upperMean90 << "-"
		<< fixed << setprecision(6) << lowerMean90 << " "
		<< fixed << setprecision(6) << upperMean95 << "-"
		<< fixed << setprecision(6) << lowerMean95 << endl;
}

vector<vector<int>> RANDOM_GRAPH(int vertices, int edges)
{
	vector<vector<int>> adjacencyMatrix;
	for (int i = 0; i < vertices; i++)
	{
		vector<int> eachVector;
		for (int j = 0; j < vertices; j++)
		{
			eachVector.push_back(0);
		}
		adjacencyMatrix.push_back(eachVector);
	}
	while (edges > 0)
	{
		int randVertex1 = rand() % vertices;
		int randVertex2 = rand() % vertices;
		if (randVertex1 != randVertex2 && adjacencyMatrix[randVertex1][randVertex2] == 0)
		{
			adjacencyMatrix[randVertex1][randVertex2] = 1;
			adjacencyMatrix[randVertex2][randVertex1] = 1;
			edges--;
		}
	}
	return adjacencyMatrix;
}

void FIXED_VERTEX_ANALYZE()
{
	vector<int> iterCounts = { 250, 500, 1000, 2500 };
	for (int iterIndex = 0; iterIndex < (int)iterCounts.size(); iterIndex++)
	{
		cout << "Statistics for ITER_COUNT=" << iterCounts[iterIndex] << endl;
		cout << "Mean     StdDev   StdErr   CL%90             CL%95" << endl;
		int vertex = 250;
		vector<int> edges = { 250, 500, 1000, 2500, 5000 };
		for (int edgeIndex = 0; edgeIndex < (int)edges.size(); edgeIndex++)
		{
			vector<float> runningTimes;
			for (int i = 0; i < iterCounts[iterIndex]; i++)
			{
				vector<Vertex> vertices;
				vector<vector<int>> graph = RANDOM_GRAPH(vertex, edges[edgeIndex]);
				int maxIndex = 0;
				int maxValue = 0;
				for (int j = 0; j < (int)graph.size(); j++)
				{
					int adj_counter = 0;
					for (int k = 0; k < (int)graph[j].size(); k++)
					{
						adj_counter += graph[j][k];
					}
					if (adj_counter > maxValue)
					{
						maxIndex = j;
						maxValue = adj_counter;
					}
					vertices.push_back(Vertex(j, adj_counter));
				}
				int root = maxIndex;

				auto startTimer = chrono::high_resolution_clock::now();
				LSC(root, vertices, graph);
				auto endTimer = chrono::high_resolution_clock::now();
				float runningTime = (endTimer - startTimer).count() / (float)1000000;
				runningTimes.push_back(runningTime);
			}
			//cout << "V=" << vertex << ", E=" << edges[edgeIndex] << endl;
			STATS_DISPLAY(runningTimes);
		}
		cout << endl;
	}
}

void FIXED_EDGE_ANALYZE()
{
	vector<int> iterCounts = { 250, 500, 1000, 2500 };
	for (int iterIndex = 0; iterIndex < (int)iterCounts.size(); iterIndex++)
	{
		cout << "Statistics for ITER_COUNT=" << iterCounts[iterIndex] << endl;
		cout << "Mean     StdDev   StdErr   CL%90             CL%95" << endl;
		vector<int> verticesCount = { 100, 200, 400, 800, 1600 };
		int edge = 1000;
		for (int vertexIndex = 0; vertexIndex < (int)verticesCount.size(); vertexIndex++)
		{
			vector<float> runningTimes;
			for (int i = 0; i < iterCounts[iterIndex]; i++)
			{
				vector<Vertex> vertices;
				vector<vector<int>> graph = RANDOM_GRAPH(verticesCount[vertexIndex], edge);
				int maxIndex = 0;
				int maxValue = 0;
				for (int j = 0; j < (int)graph.size(); j++)
				{
					int adj_counter = 0;
					for (int k = 0; k < (int)graph[j].size(); k++)
					{
						adj_counter += graph[j][k];
					}
					if (adj_counter > maxValue)
					{
						maxIndex = j;
						maxValue = adj_counter;
					}
					vertices.push_back(Vertex(j, adj_counter));
				}
				int root = maxIndex;

				auto startTimer = chrono::high_resolution_clock::now();
				LSC(root, vertices, graph);
				auto endTimer = chrono::high_resolution_clock::now();
				float runningTime = (endTimer - startTimer).count() / (float)1000000;
				runningTimes.push_back(runningTime);
			}
			STATS_DISPLAY(runningTimes);
		}
		cout << endl;
	}
}

void QUALITY_RATIO_BOUND()
{
	vector<int> iterCounts = { 10, 100, 1000, 10000 };
	for (int iterIndex = 0; iterIndex < (int)iterCounts.size(); iterIndex++)
	{
		cout << "Statistics for ITER_COUNT=" << iterCounts[iterIndex] << endl;
		vector<int> verticesCount = { 6, 8, 10, 12, 14 };
		vector<int> edgesCount = { 9, 13, 17, 21, 25 };
		vector<float> qualities;
		for (int index = 0; index < (int)verticesCount.size(); index++)
		{
			float allQualities = 0.0;
			for (int i = 0; i < iterCounts[iterIndex]; i++)
			{
				vector<Vertex> vertices;
				vector<vector<int>> graph = RANDOM_GRAPH(verticesCount[index], edgesCount[index]);
				int maxIndex = 0;
				int maxValue = 0;
				for (int j = 0; j < (int)graph.size(); j++)
				{
					int adj_counter = 0;
					for (int k = 0; k < (int)graph[j].size(); k++)
					{
						adj_counter += graph[j][k];
					}
					if (adj_counter > maxValue)
					{
						maxIndex = j;
						maxValue = adj_counter;
					}
					vertices.push_back(Vertex(j, adj_counter));
				}
				int root = maxIndex;
				LSC(root, vertices, graph);
				int finish_max = 0;
				for (int i = 1; i < (int)vertices.size(); i++)
				{
					if (vertices[i].finish > finish_max)
					{
						finish_max = vertices[i].finish;
					}
				}
				vector<Stack<int>> paths;
				Stack<int> currentStack;
				int mostDegree = 0;
				int mostTotal = 0;
				for (int i = 0; i < (int)graph.size(); i++)
				{
					int currentTotal = 0;
					for (int j = 0; j < (int)graph[0].size(); j++)
					{
						currentTotal += graph[i][j];
					}
					if (currentTotal > mostTotal)
					{
						mostDegree = i;
						mostTotal = currentTotal;
					}
				}
				int last = -1;
				vector<int> isVisited;
				for (int i = 0; i < verticesCount[index]; i++)
				{
					isVisited.push_back(0);
				}
				currentStack.makeEmpty();
				currentStack.push(mostDegree);
				LSC_BRUTEFORCE(mostDegree, last, isVisited, graph, paths, currentStack);
				int longestPathLength = 0;
				for (int i = 0; i < (int)paths.size(); i++)
				{
					if (longestPathLength < paths[i].getSize())
					{
						longestPathLength = paths[i].getSize();
					}
				}
				if (finish_max != 0 && longestPathLength != 0)
				{
					allQualities += ((float)finish_max / (float)(longestPathLength));
				}
			}
			qualities.push_back(allQualities / iterCounts[iterIndex]);
		}
		for (int i = 0; i < qualities.size(); i++)
		{
			cout << fixed << setprecision(6) << qualities[i] << endl;
		}
	}
}

void CORRECTNESS_RATIO_BOUND()
{
	int testCount = 15;
	int verticesCount = 8;
	int edgesCount = 20;
	vector<vector<vector<int>>> failed_graph;
	vector<vector<Vertex>> failed_heuristic;
	vector<Stack<int>> failed_bruteforce;
	for (int i = 0; i < testCount; i++)
	{
		vector<Vertex> vertices;
		vector<vector<int>> graph = RANDOM_GRAPH(verticesCount, edgesCount);
		int maxIndex = 0;
		int maxValue = 0;
		for (int j = 0; j < (int)graph.size(); j++)
		{
			int adj_counter = 0;
			for (int k = 0; k < (int)graph[j].size(); k++)
			{
				adj_counter += graph[j][k];
			}
			if (adj_counter > maxValue)
			{
				maxIndex = j;
				maxValue = adj_counter;
			}
			vertices.push_back(Vertex(j, adj_counter));
		}
		int root = maxIndex;
		LSC(root, vertices, graph);
		int finish_max = 0;
		for (int i = 1; i < (int)vertices.size(); i++)
		{
			if (vertices[i].finish > finish_max)
			{
				finish_max = vertices[i].finish;
			}
		}
		vector<Stack<int>> paths;
		Stack<int> currentStack;
		int mostDegree = 0;
		int mostTotal = 0;
		for (int i = 0; i < (int)graph.size(); i++)
		{
			int currentTotal = 0;
			for (int j = 0; j < (int)graph[0].size(); j++)
			{
				currentTotal += graph[i][j];
			}
			if (currentTotal > mostTotal)
			{
				mostDegree = i;
				mostTotal = currentTotal;
			}
		}
		int last = -1;
		vector<int> isVisited;
		for (int i = 0; i < verticesCount; i++)
		{
			isVisited.push_back(0);
		}
		currentStack.makeEmpty();
		currentStack.push(mostDegree);
		LSC_BRUTEFORCE(mostDegree, last, isVisited, graph, paths, currentStack);
		int longestPathLength = 0;
		Stack<int> bruteforce_path;
		for (int i = 0; i < (int)paths.size(); i++)
		{
			if (longestPathLength < paths[i].getSize())
			{
				bruteforce_path = paths[i];
				longestPathLength = paths[i].getSize();
			}
		}
		if (finish_max != longestPathLength - 1)
		{
			failed_graph.push_back(graph);
			failed_bruteforce.push_back(bruteforce_path);
			failed_heuristic.push_back(vertices);
		}
	}
	for (int i = 0; i < (int)failed_graph.size(); i++)
	{
		cout << "Adjacency Matrix:                  ";
		for (int j = 0; j < (int)failed_graph[i].size(); j++)
		{
			if (j != 0)
			{
				cout << "                                   ";
			}
			for (int k = 0; k < (int)failed_graph[i][j].size(); k++)
			{
				cout << " " << failed_graph[i][j][k];
			}
			cout << endl;
		}
		cout << "Brute Force Algorithm Longest Path:";
		while (!failed_bruteforce[i].isEmpty())
		{
			int temp = failed_bruteforce[i].pop();
			cout << " " << temp;
		}
		cout << endl;
		int maxIndex = 0;
		int maxValue = 0;
		for (int j = 0; j < (int)failed_heuristic[i].size(); j++)
		{
			if (failed_heuristic[i][j].adj_counter > maxValue)
			{
				maxIndex = j;
				maxValue = failed_heuristic[i][j].adj_counter;
			}
		}
		int root = maxIndex;
		cout << "Heuristic Algorithm Longest Path:   ";
		PRINT_CYCLE(root, failed_heuristic[i]);
		cout << endl;
	}
}

void CORNER_TEST_CASES()
{
	vector<vector<int>> test_case_1_graph =
	{
		{ 0, 0, 0, 0, 0, 0, 0, 1 },
		{ 0, 0, 1, 0, 0, 0, 0, 0 },
		{ 0, 1, 0, 1, 0, 0, 0, 0 },
		{ 0, 0, 1, 0, 1, 0, 0, 0 },
		{ 0, 0, 0, 1, 0, 1, 0, 0 },
		{ 0, 0, 0, 0, 1, 0, 1, 0 },
		{ 0, 0, 0, 0, 0, 1, 0, 1 },
		{ 1, 0, 0, 0, 0, 0, 1, 0 },
	};
	int test_case_1_root = 2;
	vector<Vertex> test_case_1_vertices = { Vertex(0, 1), Vertex(1, 1), Vertex(2, 2), Vertex(3, 2), Vertex(4, 2), Vertex(5, 2), Vertex(6, 2), Vertex(7, 2) };
	LSC(test_case_1_root, test_case_1_vertices, test_case_1_graph);
	int finish_max = 0;
	for (int i = 1; i < (int)test_case_1_vertices.size(); i++)
	{
		if (test_case_1_vertices[i].finish > finish_max)
		{
			finish_max = test_case_1_vertices[i].finish;
		}
	}
	cout << finish_max << endl;

	vector<vector<int>> test_case_2_graph =
	{
		{ 0, 0, 0, 0, 0, 0, 1, 0, 0, 1 },
		{ 0, 0, 1, 1, 1, 1, 0, 0, 0, 0 },
		{ 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 },
		{ 0, 1, 1, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 1, 0, 0, 0, 1, 0, 0, 0, 0 },
		{ 0, 1, 0, 0, 1, 0, 0, 0, 0, 0 },
		{ 1, 0, 0, 0, 0, 0, 0, 1, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 1, 0, 1, 0 },
		{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 1 },
		{ 1, 0, 0, 0, 0, 0, 0, 0, 1, 0 },
	};
	int test_case_2_root = 1;
	vector<Vertex> test_case_2_vertices = { Vertex(0, 2), Vertex(1, 4), Vertex(2, 2), Vertex(3, 2), Vertex(4, 2), Vertex(5, 2), Vertex(6, 2), Vertex(7, 2), Vertex(8, 2), Vertex(9, 2) };
	LSC(test_case_2_root, test_case_2_vertices, test_case_2_graph);
	finish_max = 0;
	for (int i = 1; i < (int)test_case_2_vertices.size(); i++)
	{
		if (test_case_2_vertices[i].finish > finish_max)
		{
			finish_max = test_case_2_vertices[i].finish;
		}
	}
	cout << finish_max << endl;

}

int main()
{
	cout << "FIXED VERTEX => V=250" << endl;
	FIXED_VERTEX_ANALYZE();
	cout << endl;
	cout << "FIXED EDGE => E=1000" << endl;
	FIXED_EDGE_ANALYZE();
	cout << endl;
	cout << "QUALITY RATIO BOUND" << endl;
	QUALITY_RATIO_BOUND();
	cout << endl;
	cout << "CORRECTNESS RATIO BOUND" << endl;
	CORRECTNESS_RATIO_BOUND();
	cout << endl;
	cout << "CORNER TEST CASES" << endl;
	CORNER_TEST_CASES();
	return 0;
}