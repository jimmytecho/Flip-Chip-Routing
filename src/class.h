#ifndef CLASS_H
#define CLASS_H
#include <vector>
#include <string>
using namespace std;

class line{
	public:
		int x1;
		int x2;
		int y1;
		int y2;	
};

class path{
 public:
	 vector <line> path;

};

class point {
	public:
		int x;
		int y;
		float index;
		float virtual_index;
};

class detail{
		public:
			int sec;
			float index;
			float virtual_index;
			float x;
			float y;
};

class Side {
	public:
		vector <point> outside;
		vector <point> inside;
		vector <point> interval;
};

class Driver {
	public:
		int x0;
		int y0;
		int xf;
		int yf;
		int index;/*original numbering*/
		
		int ring;   
		int index_in_ring;/*the number in () is the index of a ring 
						    			ex: 1   4(1)   7
						        			6   2(2)   8
						        			3   5(3)   9*/ 
	
		float virtual_index;
		bool real=true;
		vector <string> path;/*the nodes to pass*/
		vector <int> pathx;
		vector <int> pathy;
		vector <int> sec;   /*sec[0]:which section this bump/driver is in, 
							  sec[1]:default "nil", if cross-sectioned, record the counterpart's section
							  1:up  2:right 3:down 4:left*/	
		int pos_a;
		int pos_b;
		vector <point> up;
		vector <point> down;
		vector <point> left;
		vector <point> right;

		bool directlink;

		vector <point> dir;	/*every 2 int becomes 1 point to go
							x, y 
							when drawing line segment, you need to cover the point*/
		vector <int> walk;	/*every 4 int becomes 1 line segment
							output file:
							wire x0 y0 xf yf
							wire x0 y0 xf yf...
							0 1 2 3, 4 5 6 7, 8 9 10 11 
							start from driver, so use driver's walk vector
							bump's walk vector remain empty/unused */

		point start;
		point end;
		vector <point> recent;

		Side sideup;
		Side sideright;
		Side sideleft;
		Side sidedown;
};

#endif
