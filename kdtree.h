/**

  * @file   kdtree.h
  * @brief Thisis a brief description.
  * @author dongjian
  * @par   Copyright (c):
  *         All Rights Reserved
  * @date   2018:04:24 
  *  @note   mattersneeding attention
  */  
#ifndef _5383BD42_370E_4C00_A25E_AD4403E5656A
#define _5383BD42_370E_4C00_A25E_AD4403E5656A
#include <vector>
#include <queue>

namespace kdtree
{
	class Point
	{
	public:
	    double x;
		double y;
		double z;

		double Distance(const Point& p) const;
		Point(double ox, double oy, double oz);
	};

	class BBox
	{
	public:
		double xmin, ymin, zmin;
		double xmax, ymax, zmax;

		bool Contains(const Point& p) const;
		bool Contains(const BBox& bbox) const;
		bool Intersects(const BBox& bbox) const;

		BBox(double x0, double x1, double y0, double y1, double z0, double z1);

		static BBox UniverseBBox();
		double ShortestDistance(const Point& p) const;
	};

	enum Dimension
	{
		X = 0,
		Y = 1,
		Z = 2
	};

	class Node
	{
	public:
		Node* left;			//left child
		Node* right;		//right child
		int begin;			//start index [close
		int end;				//end index  (open
		Dimension dimension;	//cut dimension
		double pivot;		//cut value

		Node(int b, int e, Dimension dim, double piv);
		bool IsLeaf() const;

		BBox LeftSubTreeEnvelope(const BBox& current_entext) const;
		BBox RightSubTreeEnvelope(const BBox& current_entext) const;
	};

	class KDTree
	{
	public:
		void SetInput(const std::vector<Point>& inputs);

		void RangeSearch(const BBox& bbox, std::vector<int>& indices) const ;

		int NearestNeighbourSearch(const Point& searchpoint) const ;

		void  NearestNeighbourSearch(const Point& searchpoint, int k, std::vector<int>& indices) const;

		static int _leaf_max_size;
	private:
		const std::vector<Point>* pData;
		std::vector<int> indices;
		Node* pRoot;

	private:
		struct comparepair
		{
			bool operator()(const std::pair<int, double>&, const std::pair<int, double>&);
		};

		typedef std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, comparepair> index_distance;

		Node* DivideTree(const std::vector<Point>& inputs,int start,int end);
		Dimension SelectDimension(const std::vector<Point>& inputs, int start, int end);
		int Partition(int start, int end, Dimension dim);

		double GetAt(int index, Dimension dim);

		void _range_search(const BBox& search_bbox, Node* pNode, const BBox& current_extent,std::vector<int>& indices) const ;

		int _nearest_neighbour(const Point& searchPoint, Node* pNode, const BBox& current_extent, double& current_shorest_dis) const ;

		void _nearest_neighbour(const Point& searchPoint, int k, index_distance& ins, Node* pNode, const BBox& current_extent, double& current_shorest_dis) const ;
	};

}
#endif //_5383BD42_370E_4C00_A25E_AD4403E5656A