#include "kdtree.h"
#include <algorithm>
#include <iostream>
using namespace std;
using namespace kdtree;

Point::Point(double ox, double oy, double oz) :
x(ox),
y(oy),
z(oz)
{};

double Point::Distance(const Point& p) const
{
	double dx = p.x - x;
	double  dy = p.y - y;
	double dz = p.z - z;
	return sqrt(dx*dx + dy*dy + dz*dz);
}

BBox Node::LeftSubTreeEnvelope(const BBox& current_entext) const
{
	BBox leftRegion(current_entext);
	switch (dimension)
	{
	case X:
		leftRegion.xmax = pivot;
		break;
	case Y:
		leftRegion.ymax = pivot;
		break;
	case Z:
		leftRegion.zmax = pivot;
		break;
	}
	return leftRegion;
}

BBox Node::RightSubTreeEnvelope(const BBox& current_entext) const
{
	BBox rightRegion(current_entext);
	switch (dimension)
	{
	case X:
		rightRegion.xmin = pivot;
		break;
	case Y:
		rightRegion.ymin = pivot;
		break;
	case Z:
		rightRegion.zmin = pivot;
		break;
	}
	return rightRegion;
}

bool BBox::Contains(const BBox& bbox) const
{
	if (bbox.xmin < xmin)		return false;
	if (bbox.xmax > xmax)		return false;
	if (bbox.ymin < ymin)		return false;
	if (bbox.ymax < ymax)		return false;
	if (bbox.zmin < zmin)		return false;
	if (bbox.zmax > zmax)		return false;
	return true;
}

bool BBox::Contains(const Point& p) const
{
	return p.x >= xmin && p.x <= xmax && p.y >= ymin && p.y <= ymax
		&& p.z >= zmin && p.z <= zmax;
}

bool BBox::Intersects(const BBox& bbox) const
{
	if (bbox.xmin > xmax || bbox.xmax < xmin) return false;
	if (bbox.ymin > xmax || bbox.ymax < ymin) return false;
	if (bbox.zmin > zmax || bbox.zmax < zmin) return false;
	return true;
}

BBox BBox::UniverseBBox()
{
	double DOUBLE_MAX = std::numeric_limits<double>::max();
	return BBox(-DOUBLE_MAX, DOUBLE_MAX, -DOUBLE_MAX, 
		DOUBLE_MAX, -DOUBLE_MAX, DOUBLE_MAX);
}

BBox::BBox(double x0, double x1, double y0, double y1, double z0, double z1) :
xmin(x0),
xmax(x1),
ymin(y0),
ymax(y1),
zmin(z0),
zmax(z1)
{}

double BBox::ShortestDistance(const Point& p) const
{
	double dx = xmin - p.x > 0 ? xmin - p.x : (p.x - xmax > 0 ? p.x - xmax : 0);
	double dy = ymin - p.y > 0 ? ymin - p.y : (p.y - ymax > 0 ? p.y - ymax : 0);
	double dz = zmin - p.z > 0 ? zmin - p.z : (p.z - zmax > 0 ? p.z - zmax : 0);
	return sqrt(dx*dx + dy*dy + dz*dz);
}


Node::Node(int b, int e, Dimension dim, double piv) :
begin(b),
end(e),
dimension(dim),
pivot(piv),
left(nullptr),
right(nullptr)
{

}

bool Node::IsLeaf() const
{
	return left == nullptr && right == nullptr;
}

Dimension KDTree::SelectDimension(const std::vector<Point>& inputs, int start, int end)
{
	struct cmpobj
	{
		Dimension dim;
		cmpobj(Dimension d) :
			dim(d){}

		bool operator()(const Point& p0, const Point p1) const
		{
			if (dim == X)
				return p0.x < p1.x;
			if (dim == Y)
				return p0.y < p1.y;
			if (dim == Z)
				return p0.z < p1.z;
			return false;
		}
	};

	double span[3];

	auto pair =  std::minmax_element(inputs.begin() + start, inputs.begin() + end,cmpobj(X));
	span[0] = pair.second->x - pair.first->x;

	pair = std::minmax_element(inputs.begin() + start, inputs.begin() + end, cmpobj(Y));
	span[1] = pair.second->y - pair.first->y;

	pair = std::minmax_element(inputs.begin() + start, inputs.begin() + end, cmpobj(Z));
	span[2] = pair.second->z - pair.first->z;

	auto index = std::distance(span, std::max_element(span, span + 3));

	return static_cast<Dimension>(index);
}

double KDTree::GetAt(int index, Dimension dim)
{
	auto p = (*pData)[indices[index]];
	return  dim == X ? p.x : (dim == Y ? p.y : p.z);
}

int KDTree::Partition(int start, int end, Dimension dimension)
{
	int size = end - start;
	if (size <= 0)
	{
		cout << "a serious error occurs " << start << "\t" << endl;
		return -1;
	}
		
	struct cmpobj
	{
		Dimension dim;
		const std::vector<Point>* pData;
		bool operator()(int i, int j) const 
		{
			if (X == dim)
				return (*pData)[i].x < (*pData)[j].x;
			if (Y == dim)
				return (*pData)[i].y < (*pData)[j].y;
			if (Z == dim)
				return (*pData)[i].z < (*pData)[j].z;

			return true;
		}
		cmpobj(Dimension dimension, const std::vector<Point>* pInputData):
			dim(dimension),
			pData(pInputData)
		{}
	};
	
	int median = start +size / 2;
	std::nth_element(indices.begin() + start, indices.begin()+median , indices.begin() + end, cmpobj(dimension, pData));
	return median;
}

int KDTree::NearestNeighbourSearch(const Point& searchpoint) const 
{
	double shortestDistance = std::numeric_limits<double>::max();
    return _nearest_neighbour(searchpoint, pRoot, BBox::UniverseBBox(), shortestDistance);
}

void KDTree::NearestNeighbourSearch(const Point& searchpoint, int k, std::vector<int>& indices) const 
{
	double shortestDistance = std::numeric_limits<double>::max();
	indices.clear();
	 
	index_distance neighbours;
	_nearest_neighbour(searchpoint, k, neighbours, pRoot, BBox::UniverseBBox(), shortestDistance);

	while (!neighbours.empty())
	{
		indices.push_back(neighbours.top().first);
		neighbours.pop();
	}
}

bool KDTree::comparepair::operator()(const std::pair<int, double>& p0, const std::pair<int, double>& p1)
{
	return p0.second < p1.second;
}

void KDTree::_nearest_neighbour(const Point& searchPoint, int k, index_distance& ins,
	Node* pNode, const BBox& current_extent, double& current_shorest_dis) const 
{
	double min_shortest_distance = current_extent.ShortestDistance(searchPoint);
	if (min_shortest_distance >= current_shorest_dis)
		return;

	if (pNode->IsLeaf())
	{
		for (auto i = pNode->begin; i < pNode->end; ++i)
		{
			double distance = (*pData)[indices[i]].Distance(searchPoint);
			if (ins.size() < k)
			{
				ins.push(pair<int, double>(indices[i], distance));	//add element
				if (k == ins.size())
					current_shorest_dis = ins.top().second;
			}
			else 
			{
				if (distance < current_shorest_dis)
				{
					ins.pop();
					ins.push(pair<int, double>(indices[i], distance));//add element
					current_shorest_dis = ins.top().second;
				}
			}
		}
	}
	else
	{
		BBox leftRegion = pNode->LeftSubTreeEnvelope(current_extent);
		BBox rightRegion = pNode->RightSubTreeEnvelope(current_extent);

		double dis_to_left = leftRegion.ShortestDistance(searchPoint);
		double dis_to_right = rightRegion.ShortestDistance(searchPoint);

		if (dis_to_left < dis_to_right)
		{
			_nearest_neighbour(searchPoint,k,ins,pNode->left,leftRegion,current_shorest_dis);
			 _nearest_neighbour(searchPoint, k,ins,pNode->right, rightRegion, current_shorest_dis);
		}
		else
		{
			_nearest_neighbour(searchPoint,k,ins, pNode->right, rightRegion, current_shorest_dis);
			_nearest_neighbour(searchPoint,k,ins, pNode->left, leftRegion, current_shorest_dis);
		}
	}
}


int KDTree::_nearest_neighbour(const Point& searchPoint, Node* pNode, const BBox& current_extent, double& current_shorest_dis) const 
{
	double min_shortest_distance = current_extent.ShortestDistance(searchPoint);
	if (min_shortest_distance >= current_shorest_dis)
		return -1;

	if (pNode->IsLeaf())
	{
		int shortestindex =-1;
		for (auto i = pNode->begin; i < pNode->end;++i)
		{
			double distance = (*pData)[indices[i]].Distance(searchPoint);
			if (distance < current_shorest_dis)
			{
				shortestindex = indices[i];
				current_shorest_dis = distance;
			}
		}
		return shortestindex;
	}
	else
	{
		BBox leftRegion = pNode->LeftSubTreeEnvelope(current_extent);
		BBox rightRegion = pNode->RightSubTreeEnvelope(current_extent);

		double dis_to_left = leftRegion.ShortestDistance(searchPoint);
		double dis_to_right = rightRegion.ShortestDistance(searchPoint);

		if (dis_to_left < dis_to_right)
		{
			int left = _nearest_neighbour(searchPoint, pNode->left, leftRegion, current_shorest_dis);
			int right = _nearest_neighbour(searchPoint, pNode->right, rightRegion, current_shorest_dis);

			return right == -1 ? left : right;
		}
		else
		{
			int right = _nearest_neighbour(searchPoint, pNode->right, rightRegion, current_shorest_dis);
			int left = _nearest_neighbour(searchPoint, pNode->left, leftRegion, current_shorest_dis);
			return left == -1 ? right : left;
		}

		return -1;
	}
}

int KDTree::_leaf_max_size = 15;

void KDTree::SetInput(const std::vector<Point>& inputs)
{
	pData = &inputs;
	indices.resize(inputs.size());
	for (int i = 0, n = inputs.size(); i < n; ++i)
		indices[i] = i;
	pRoot = DivideTree(inputs, 0, inputs.size());
}

Node* KDTree::DivideTree(const std::vector<Point>& inputs, int start, int end)
{
	//cout << "build " << start << "\t" << end << endl;

	int size = end - start;
	if (size <= 0)
		return nullptr;

	Dimension dim = SelectDimension(inputs, start, end);
	int median = Partition(start, end, dim);

	Node* pNode = new Node(start, end, dim, GetAt(median, dim));

	if (size > _leaf_max_size)
	{
		pNode->left = DivideTree(inputs, start, median);
		pNode->right = DivideTree(inputs, median, end);
	}

	return pNode;
}

void KDTree::RangeSearch(const BBox& bbox, std::vector<int>& indices) const 
{
	BBox universe_bbox = BBox::UniverseBBox();
	_range_search(bbox, pRoot, universe_bbox, indices);
}


void KDTree::_range_search(const BBox& search_bbox, Node* pNode, const BBox& current_extent, std::vector<int>& ins) const 
{
	if (nullptr == pNode)
		return;

	if (pNode->IsLeaf())
	{
		for (int i = pNode->begin; i < pNode->end; ++i)
		{
			const Point& p = (*pData)[indices[i]];
			if (search_bbox.Contains(p))
				ins.push_back(indices[i]);
		}
	}
	else
	{
		//trim bounding box
		BBox leftRegion = pNode->LeftSubTreeEnvelope(current_extent);

		if (search_bbox.Contains(leftRegion))
		{
			for (int i = pNode->left->begin; i < pNode->left->end; ++i)
				ins.push_back(indices[i]);
		}
		else if (search_bbox.Intersects(leftRegion))
		{
			_range_search(search_bbox, pNode->left, leftRegion, ins);
		}

		BBox rightRegion = pNode->RightSubTreeEnvelope(current_extent);

		if (search_bbox.Contains(rightRegion))
		{
			for (int i = pNode->right->begin; i < pNode->right->end; ++i)
				ins.push_back(indices[i]);
		}
		else if (search_bbox.Intersects(rightRegion))
		{
			_range_search(search_bbox, pNode->right, rightRegion, ins);
		}
	}
}
