/*
 * vol2mesh.c
 *
 * Author: A. Berariu 
 * Date: Oct 08, 2009
 * Update: Oct 22, 2014 
 */

#include <stdio.h>
#include <stdlib.h>

/**
 * @brief This function is used to convert the volume-based topology (of Fire input files) 
 * to unstructured mesh topology and geometry (for vtk format output).
 *
 * The volume topology is received via the input array LCC (link cell-cell). One LCC[nghbId][cellId]
 * entry in the array stores the neighboring cell of cellId in direction nghbId.
 * Each volume cell is considered to be a 1x1x1 cube. The transformation result is stored in three 
 * variables:
 *      elems: 2D array storing the mesh elements. An element is defined by 8 points(nodes):
 *             elems[i][NC] -> i-th node of mesh element NC (see also scheme below)
 *      points: 2D array storing the (x,y,z) coordinates for each point (node):
 *              points[0][i] -> x coord of point i; points[1][i]->y coord of point i; ...
 *      nodeCnt: the total number of points in the mesh (size of second dimension of points array)
 *
 * @param startIdx: the index of the first entry in LCC
 * @param endIdx: the index of the last entry in LCC
 * @param LCC: the volume cells topology
 *
 * return:
 *      @param nodeCnt: the total number of points (after conversion)
 *      @param points: the array containing the coordinate of the points
 *      @param elems: the array containing the mesh elements
 */
int vol2mesh(int startIdx, int endIdx, int** LCC, int* nodeCnt, int*** points,	int*** elems)
{

	int elemCnt = endIdx - startIdx + 1;

	/************************************************
	 * Algorithm for volume elements to mesh elements conversion
	 *
	 *  - vol. elements are defined by their 6 neighbors
	 *  - (hexahedral) mesh elements are defined by 8 nodes
	 *  - each node is shared by 4 mesh elements
	 *
	 *
	 *     for each volume element (in their address order) do
	 *          for each node n (out of 8) do
	 *               - check if one of the neighbors sharing this node
	 *                    has a smaller address (aka has already been processed)
	 *               - if so, then copy the node definition
	 *               - otherwise, define this new node
	 * */

	/** map from current node idx to the node idxs in each of the 3 neighbors
	 *                                  ^ z
	 *           5--------6             |
	 *          /|       /|             |
	 *         4-|------7 |             |
	 *         | |      | |             ---------->
	 *         | 1------|-2            /           y
	 *         |/       |/            /
	 *         0--------3            /
	 *                              `´x
	 *           H N
	 *           |/
	 *       W - * - E
	 *          /|
	 *         S L
	 *
	 *  neighboring direction:
	 *              S - 0
	 *              E - 1
	 *              N - 2
	 *              W - 3
	 *              L - 4
	 *              H - 5
	 * */
	int node2nghbMap[8][6] = { {0, 1, 3, 3, 4, 4}, // node 7 on S nghb, node 4 on W nghb, node 5 on L nghb
			{3, 2, 4, 5, 2, 0}, {4, 6, 1, 1, 2, 3}, {0, 2, 4, 7, 1, 0}, {3, 7,
					0, 5, 5, 0}, {3, 6, 5, 1, 2, 4}, {5, 2, 2, 7, 1, 5}, {0, 6,
					5, 3, 1, 4}};

	/**
	 * define the point coordinates along with the cell connectivity
	 *
	 * for each newly defined node, set the coordinates related to the
	 * previous node in the cell:
	 *
	 * node 0: ********** ?
	 * node 1: x1 = x0 - 1
	 * node 2: y2 = y1 + 1
	 * node 3: x3 = x2 + 1
	 * node 4: z4 = z3 + 1
	 *         y4 = y3 - 1
	 * node 5: x5 = x4 - 1
	 * node 6: y6 = y5 + 1
	 * node 7: x7 = x6 + 1
	 *
	 * in order to find the coordinates of a node with index 0,
	 * search a nghb with smaller cell index and derive the coordinates from one of its points:
	 *      search among the nghbs which do not share this node (otherwise, it would have
	 * already be defined):
	 *              E: y0 = y0E - 1
	 *              H: z0 = z0H - 1
	 *              N: x0 = x0N + 1
	 * */
	//coordinate offsets in own element/cell
	int OFFSET[8][3] = { {0, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {1, 0, 0},
										{0, -1, 1}, {-1, 0, 0}, {0, 1, 0}, {1, 0, 0}};
	// the points (coordinate):
	//      [0] -> x coord
	//      [1] -> y coord
	//      [2] -> z coord
	if ((*points = (int **) malloc(3 * sizeof(int*))) == NULL)
	{
		printf("malloc() POINTS 1st dim. failed\n");
		return -1;
	}

	int i;
	for (i = 0; i < 3; i++)
	{ //note: oversized array .. better safe than sorry ;)
		if (((*points)[i] = (int *) malloc((elemCnt + 1) * sizeof(int) * 8))
				== NULL)
		{
			printf("malloc() POINTS 2nd dim. failed\n");
			return -1;
		}
	}

	// the cells (with node indecies)
	//int ***elems; // LCC's dual
	if (((*elems) = (int **) malloc(8 * sizeof(int*))) == NULL)
	{
		printf("malloc() ELEMS 1st dim. failed\n");
		return -1;
	}
	for (i = 0; i < 8; i++)
	{
		if (((*elems)[i] = (int *) malloc((elemCnt + 1) * sizeof(int))) == NULL)
		{
			printf("malloc() ELEMS 2nd dim. failed\n");
			return -1;
		}
	}

	int nodes = 0;
	int nodeIdx = 0;
	int nghbIdx = 0;
	int coordIdx = 0;
	int currNghb = 0;
	int foundVal = 0; // aka false

	/*******Start iteration loop****************************************/
	int NC;
	for (NC = startIdx; NC <= endIdx; NC++)
	{ // for each volume cell

		for (nodeIdx = 0; nodeIdx < 8; nodeIdx++)
		{ // for each node
			foundVal = 0;
			//check the 3 neighbours
			for (nghbIdx = 0; nghbIdx < 6; nghbIdx = nghbIdx + 2)
			{
				currNghb = LCC[NC][node2nghbMap[nodeIdx][nghbIdx]];
				if (currNghb < NC)
				{
					(*elems)[nodeIdx][NC] =
							(*elems)[node2nghbMap[nodeIdx][nghbIdx + 1]][currNghb];
					foundVal = 1; //aka true
					break;
				}
			}
			// if the node is not defined yet, give it a new address (id)
			if (!foundVal)
			{
				nodes++;
				(*elems)[nodeIdx][NC] = nodes;
				// set the point coordinates
				if (nodes == 1)
				{ // the first node ever :)
					(*points)[0][nodes] = (*points)[1][nodes] = (*points)[2][nodes] = 0;
				} else if (nodeIdx == 0)
				{ //special case, search for coordinates by neighbors
					// E nghb
					currNghb = LCC[NC][1];
					if (NC > currNghb)
					{
						(*points)[0][nodes] = (*points)[0][currNghb];
						(*points)[1][nodes] = (*points)[1][currNghb] - 1;
						(*points)[2][nodes] = (*points)[2][currNghb];
					}
					else
					{
						currNghb = LCC[NC][2];
						if (NC > currNghb)
						{
							(*points)[0][nodes] = (*points)[0][currNghb] + 1;
							(*points)[1][nodes] = (*points)[1][currNghb];
							(*points)[2][nodes] = (*points)[2][currNghb];
						}
						else
						{
							currNghb = LCC[NC][5];
							if (NC > currNghb)
							{
								(*points)[0][nodes] = (*points)[0][currNghb];
								(*points)[1][nodes] = (*points)[1][currNghb];
								(*points)[2][nodes] = (*points)[2][currNghb] - 1;
							}
						}
					}
				}
				else
				{ // all other node coord. are determined from own node coord
					for (coordIdx = 0; coordIdx < 3; coordIdx++)
						(*points)[coordIdx][nodes] = (*points)[coordIdx][(*elems)[nodeIdx - 1][NC]] + OFFSET[nodeIdx][coordIdx];
				}
			}
		} // for each node

	} // end for each volume cell/mesh element

	*nodeCnt = nodes;

	return 0;
}
