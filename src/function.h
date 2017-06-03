#ifndef FUNCTION_H
#define FUNCTION_H

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "class.h"
#include <cmath>

void link(vector < vector <Driver> > &graph, vector <path> &mypath, int i, int j, int width, int space) {
	int k, n, id, sp;
	int a = 0;
	int skip = -1;
	int skip2 = -1;
	bool sskip = true;
	sp = width / 2 + space;
	vector <point> A, B;// A1, B1;
	B = graph[i][j].sideup.interval;
	A = graph[i][j].sidedown.interval;
	//set n;
	n = B.size();
	if (n > A.size()) {
		n = A.size();
	}
	//check if there is wierd points
	if (n != 0) {
		//sorting
		for (int k3 = 0; k3 < n; k3++) {
			point temp = A[k3];
			int k4 = k3 - 1;
			if (k4 > -1) {
				while (k4 > -1 && (temp.x > A[k4].x)) {
					A[k4 + 1] = A[k4];
					A[k4] = temp;
					k4--;
				}
			}
		}
		for (int k3 = 0; k3 < n; k3++) {
			point temp = B[k3];
			int k4 = k3 - 1;
			if (k4 > -1) {
				while (k4>-1 && (temp.x > B[k4].x)) {
					B[k4 + 1] = B[k4];
					B[k4] = temp;
					k4--;
				}
			}
		}
		//do from right!!!!!!
		vector <point> temp;
		vector <point> prev;
		vector <point>::iterator itP, itP2;
		//link A to B;
		for (k = 0; k < n; k++) {
			//the first line
			if (k == 0) {
				//to right
				if (B[k].x > A[k].x) {
					point T;
					T.x = A[k].x;
					T.y = A[k].y;
					temp.push_back(T);
					T.x = B[k].x;
					temp.push_back(T);
					T.y = B[k].y;
					temp.push_back(T);
				}
				//to left or up
				else {
					point T;
					T.x = A[k].x;
					T.y = A[k].y;
					temp.push_back(T);
					T.y = B[k].y;
					temp.push_back(T);
					//to left
					if (B[k].x < A[k].x) {
						T.x = B[k].x;
						temp.push_back(T);
					}
				}
			}
			// the second to last line;
			else {
				//to right		
				if (B[k].x > A[k].x) {
					point T;
					T.x = A[k].x;
					T.y = A[k].y;
					temp.push_back(T);
					//start crawling;
					for (itP = prev.begin(); itP < prev.end(); itP++) {
						if ((*itP).x - sp < B[k].x) {
							if ((*itP).x - sp > T.x) {
								T.x = (*itP).x - sp;
								temp.push_back(T);
							}
							if ((*itP).y + sp > T.y) {
								T.y = (*itP).y + sp;
								temp.push_back(T);
							}
						}
					}
					//empty prev after used;
					prev.clear();
					if (T.x != B[k].x) {
						T.x = B[k].x;
						temp.push_back(T);
					}
					T.y = B[k].y;
					temp.push_back(T);
				}
				//to left or to up;
				else {
					point T;
					T.x = A[k].x;
					T.y = A[k].y;
					temp.push_back(T);
					int counter = 0;
					//start crawling; //if to up no crawling will occur;
					for (itP = prev.begin(); itP < prev.end(); itP++) {
						if ((*itP).x - sp < T.x) {
							if (counter == 0) {
								T.y = (*itP).y - sp;
								temp.push_back(T);
								T.x = (*itP).x - sp;
								temp.push_back(T);
								counter = 1;
							}
							else {
								if (T.y < (*itP).y - sp) {
									T.y = (*itP).y - sp;
									temp.push_back(T);
								}
								if (T.x > (*itP).x - sp){
									T.x = (*itP).x - sp;
									temp.push_back(T);
								}
							}
						}
					}
					T.y = B[k].y;
					temp.push_back(T);
					//to left;
					if (T.x > B[k].x) {
						T.x = B[k].x;
						temp.push_back(T);
					}
					//empty prev after used;
					prev.clear();
				}
			}
			//record temp to prev and temp to mypath;
			for (itP = temp.begin(); itP < temp.end(); itP++) {
				prev.push_back(*itP);
			}
			//record temp to mypath;
			id = A[k].index;
			for (int k5 = 0; k5<temp.size() - 1; k5++) {
				line L;
				L.x1 = temp[k5].x;
				L.y1 = temp[k5].y;
				L.x2 = temp[k5 + 1].x;
				L.y2 = temp[k5 + 1].y;
				mypath[id].path.push_back(L);
			}
			//clear temp after use;
			temp.clear();
		}
	}
}

void link2(vector < vector <Driver> > &graph, vector <path> &mypath, int i, int j, int width, int space) {
	int k, n, id, sp;
	int a = 0;
	int skip = -1;
	int skip2 = -1;
	bool sskip = true;
	sp = width / 2 + space;
	vector <point> A, B;// , A1, B1;
	B = graph[i][j].sideleft.interval;
	A = graph[i][j].sideright.interval;
	//set n;
	n = B.size();
	if (n > A.size()) {
		n = A.size();
	}
	if (n != 0) {
		//sorting
		for (int k3 = 0; k3 < n; k3++) {
			point temp = A[k3];
			int k4 = k3 - 1;
			if (k4 > -1) {
				while (k4 > -1 && (temp.y > A[k4].y)) {
					A[k4 + 1] = A[k4];
					A[k4] = temp;
					k4--;
				}
			}
		}
		for (int k3 = 0; k3 < n; k3++) {
			point temp = B[k3];
			int k4 = k3 - 1;
			if (k4 > -1) {
				while (k4>-1 && (temp.y > B[k4].y)) {
					B[k4 + 1] = B[k4];
					B[k4] = temp;
					k4--;
				}
			}

		}
		// A right B left link A to B as if bottom to top
		vector <point> temp;
		vector <point> prev;
		vector <point>::iterator itP, itP2;
		//link A to B;
		for (k = 0; k < n; k++) {
			//the first line
			if (k == 0) {
				//to right
				if (B[k].y > A[k].y) {
					point T;
					T.x = A[k].x;
					T.y = A[k].y;
					temp.push_back(T);
					T.y = B[k].y;
					temp.push_back(T);
					T.x = B[k].x;
					temp.push_back(T);
				}
				//to left or up
				else {
					point T;
					T.x = A[k].x;
					T.y = A[k].y;
					temp.push_back(T);
					T.x = B[k].x;
					temp.push_back(T);
					//to left
					if (B[k].y < A[k].y) {
						T.y = B[k].y;
						temp.push_back(T);
					}
				}
			}
			// the second to last line;
			else {
				//to right		
				if (B[k].y > A[k].y) {
					point T;
					T.x = A[k].x;
					T.y = A[k].y;
					temp.push_back(T);
					//start crawling;
					for (itP = prev.begin(); itP < prev.end(); itP++) {
						if ((*itP).y - sp < B[k].y) {
							if (T.y < (*itP).y - sp) {
								T.y = (*itP).y - sp;
								temp.push_back(T);
							}
							if (T.x >(*itP).x - sp){
								T.x = (*itP).x - sp;
								temp.push_back(T);
							}
						}
					}
					//empty prev after used;
					prev.clear();
					if (T.y != B[k].y) {
						T.y = B[k].y;
						temp.push_back(T);
					}
					T.x = B[k].x;
					temp.push_back(T);
				}
				//to left or to up;
				else {
					point T;
					T.x = A[k].x;
					T.y = A[k].y;
					temp.push_back(T);
					int counter = 0;
					//start crawling; //if to up no crawling will occur;
					for (itP = prev.begin(); itP < prev.end(); itP++) {
						if ((*itP).y - sp < T.y) {
							if (counter == 0) {
								T.x = (*itP).x + sp;
								temp.push_back(T);
								T.y = (*itP).y - sp;
								temp.push_back(T);
								counter = 1;
							}
							else {
								if (T.x >(*itP).x + sp) {
									T.x = (*itP).x + sp;
									temp.push_back(T);
								}
								if (T.y > (*itP).y - sp) {
									T.y = (*itP).y - sp;
									temp.push_back(T);
								}
							}
						}
					}
					T.x = B[k].x;
					temp.push_back(T);
					//to left;
					if (T.y > B[k].y) {
						T.y = B[k].y;
						temp.push_back(T);
					}
					//empty prev after used;
					prev.clear();
				}
			}
			//record temp to prev and temp to mypath;
			for (itP = temp.begin(); itP < temp.end(); itP++) {
				prev.push_back(*itP);
			}
			//record temp to mypath;
			id = A[k].index;
			for (int k5 = 0; k5<temp.size() - 1; k5++) {
				line L;
				L.x1 = temp[k5].x;
				L.y1 = temp[k5].y;
				L.x2 = temp[k5 + 1].x;
				L.y2 = temp[k5 + 1].y;
				mypath[id].path.push_back(L);
			}
			//clear temp after use;
			temp.clear();
		}
	}
}

void initialize(vector<Driver> &driver, vector<Driver> &bump, vector<int> &indexd, vector<int> &indexb, int sec, int wi, int le) {
	int i, x, y, z, z1;
	for (i = 0; i < driver.size(); i++) {
		//first point;
		x = (driver[i].xf + driver[i].x0) / 2;
		y = (driver[i].yf + driver[i].y0) / 2;
		driver[i].start.x = x;
		driver[i].start.y = y;
	}
	//check if it is direct link;
	if (sec == 1) {
		for (i = 0; i < driver.size(); i++) {
			z = driver[i].index;
			if (bump[indexb[z]].pos_a == 1) {
				bump[indexb[z]].directlink = true;
				if (i != driver.size() - 1) {
					z1 = driver[i + 1].index;
					if (bump[indexb[z1]].pos_a == 1) {
						if (bump[indexb[z1]].x0 < bump[indexb[z]].x0) bump[indexb[z]].directlink = false;
					}
				}
			}
			else bump[indexb[z]].directlink = false;
		}
	}
	else if (sec == 2) {
		for (i = 0; i < driver.size(); i++) {
			z = driver[i].index;
			if (bump[indexb[z]].pos_b == wi - 2) {
				bump[indexb[z]].directlink = true;
				if (i != driver.size()) {
					z1 = driver[i + 1].index;
					if (bump[indexb[z1]].pos_b == wi - 2) {
						if (bump[indexb[z1]].y0 > bump[indexb[z]].y0) bump[indexb[z]].directlink = false;
					}
				}
			}
			else bump[indexb[z]].directlink = false;
		}
	}
	else if (sec == 3) {
		for (i = 0; i < driver.size(); i++) {
			z = driver[i].index;
			if (bump[indexb[z]].pos_a == le - 2) {
				bump[indexb[z]].directlink = true;
				if (i != driver.size()) {
					z1 = driver[i + 1].index;
					if (bump[indexb[z1]].pos_a == le - 2) {
						if (bump[indexb[z1]].x0 > bump[indexb[z]].x0) bump[indexb[z]].directlink = false;
					}
				}
			}
			else bump[indexb[z]].directlink = false;
		}
	}
	else if (sec == 4) {
		for (i = 0; i < driver.size(); i++) {
			z = driver[i].index;
			if (bump[indexb[z]].pos_b == 1) {
				bump[indexb[z]].directlink = true;
				if (i != driver.size()) {
					z1 = driver[i + 1].index;
					if (bump[indexb[z1]].pos_b == 1) {
						if (bump[indexb[z1]].y0 < bump[indexb[z]].y0) bump[indexb[z]].directlink = false;
					}
				}
			}
			else bump[indexb[z]].directlink = false;
		}
	}
}

void Walk(vector < vector <Driver> > &graph, vector<Driver> &driver, vector<Driver> &bump, vector<int> &indexd, vector<int> &indexb, vector <path> &mypath, int width, int space) {

	int sp = width / 2 + space;//spacing = certain value;
	int a, i, j, k, t, ii, z;
	int position;
	vector<Driver>::iterator itD;

	//bump[indexb[i]].path.end() == bump[indexb[i]].path.begin() change to 
	//bump[indexb[i]].pos_a == 1
	//determine end;
	i = 0;
	a = -1;
	while (i < driver.size()) {
		//if this is a direct link
		a = driver[i].index;//a is the index
		if ((bump[indexb[a]].pos_a == 1) && (bump[indexb[a]].directlink)) {
			driver[i].end.x = (bump[indexb[a]].x0 + bump[indexb[a]].xf) / 2;
			driver[i].end.y = (bump[indexb[a]].y0 + bump[indexb[a]].yf) / 2;
			i++;
		}
		else {
			j = i;
			t = 0;
			bool condition = true;
			while (condition) {
				j++;
				t++;
				z = driver[j].index;
				if ((bump[indexb[z]].pos_a == 1) && (bump[indexb[z]].directlink)) {
					condition = false;
				}
				else if (j == driver.size()) {
					condition = false;
				}
				//record the other drivers that goes to the same place as i , there will be t of them, and the next is j(j is not one of them)
			}
			position = bump[indexb[z]].pos_b;
			//back from j
			if (j != driver.size()) {
				for (k = 0; k < t; k++) {
					driver[j - k - 1].end.x = bump[indexb[z]].x0 - (k + 1)* sp;//incase suitcase does not work
					vector <point> suitcase;
					suitcase = graph[1][position - 1].sidedown.interval;
					for (int su = 0; su < suitcase.size(); su++) {
						if (suitcase[su].index == driver[j - k - 1].index) driver[j - k - 1].end.x = suitcase[su].x;
					}
					driver[j - k - 1].end.y = bump[indexb[z]].yf;
					driver[j - k - 1].end.index = driver[j - k - 1].index;
					graph[1][position - 1].sideup.interval.push_back(driver[j - k - 1].end);
			}
				i = j;
			}
			else {
				a = driver[i - 1].index;//need to use the previous one;
				position = bump[indexb[a]].pos_b;
				for (k = 0; k < t; k++) {
					driver[i + k].end.x = bump[indexb[a]].xf + (k + 1)* sp;
					vector <point> suitcase;
					suitcase = graph[1][position + 1].sidedown.interval;
					for (int su = 0; su < suitcase.size(); su++) {
						if (suitcase[su].index == driver[i + k].index) driver[i + k].end.x = suitcase[su].x;
					}
					driver[i + k].end.y = bump[indexb[a]].yf;
					driver[i + k].end.index = driver[i + k].index;
					graph[1][position + 1].sideup.interval.push_back(driver[i + k].end);
				}
				i = j;
			}
		}
	}
	//link start to end
	k = 0;
	a = 1;
	for (ii = 0; ii < driver.size(); ii++) {
		if (a == 1) {
			if ((driver[ii].start.x >= driver[ii].end.x) || (ii + 1 == driver.size())) {
				//i -1 back to k;	//draw i-1
				//create recent;
				if (ii + 1 != driver.size()) i = ii - 1;
				else i = ii;

				point A;
				A.x = driver[i].start.x;
				A.y = driver[i].y0 - sp;
				driver[i].recent.push_back(A);
				A.x = driver[i].end.x;
				driver[i].recent.push_back(A);
				A.y = driver[i].end.y;
				driver[i].recent.push_back(A);
				//record to dir		
				A = driver[i].start;
				driver[i].dir.push_back(A);
				vector <point>::iterator itP;
				for (itP = driver[i].recent.begin(); itP < driver[i].recent.end(); itP++) {
					driver[i].dir.push_back(*itP);
				}
				// i-1 to k
				for (t = i - 1; t >= k; t--) {
					A.x = driver[t].start.x;
					A.y = driver[t].y0 - sp;
					driver[t].recent.push_back(A);
					if (driver[t + 1].start.x - sp < driver[t].end.x) {
						A.x = driver[t + 1].start.x - sp;
						driver[t].recent.push_back(A);
					}
					for (itP = driver[t + 1].recent.begin(); itP < driver[t + 1].recent.end(); itP++) {

						if ((*itP).x - sp < driver[t].end.x) {
							A.x = (*itP).x - sp;
							A.y = (*itP).y - sp;
							driver[t].recent.push_back(A);
						}
					}
					if (A.x != driver[t].end.x) {
						A.x = driver[t].end.x;
						driver[t].recent.push_back(A);
					}
					A.y = driver[t].end.y;
					driver[t].recent.push_back(A);
					//record to dir
					A = driver[t].start;
					driver[t].dir.push_back(A);
					vector <point>::iterator itP;
					for (itP = driver[t].recent.begin();
					itP < driver[t].recent.end(); itP++) {
						driver[t].dir.push_back(*itP);
					}
				}
				k = i + 1;
				a = 0;
			}
		}
		else if (a == 0) {
			if ((driver[ii].start.x <= driver[ii].end.x) || (ii + 1 == driver.size())) {
				//k to i;
				//draw k;
				if (ii + 1 != driver.size()) i = ii - 1;
				else i = ii;
				point A;
				A.x = driver[k].start.x;
				A.y = driver[k].y0 - sp;
				driver[k].recent.push_back(A);
				A.x = driver[k].end.x;
				driver[k].recent.push_back(A);
				A.y = driver[k].end.y;
				driver[k].recent.push_back(A);
				//record to dir
				A = driver[k].start;
				driver[k].dir.push_back(A);

				vector <point>::iterator itP;
				for (itP = driver[k].recent.begin(); itP < driver[k].recent.end(); itP++) {
					driver[k].dir.push_back(*itP);
				}
				//t = k+1 to i
				for (t = k + 1; t <= i; t++) {
					A.x = driver[t].start.x;
					A.y = driver[t].y0 - sp;
					driver[t].recent.push_back(A);
					if (driver[t - 1].start.x + sp > driver[t].end.x) {
						A.x = driver[t - 1].start.x + sp;
						driver[t].recent.push_back(A);
					}
					for (itP = driver[t - 1].recent.begin(); itP < driver[t - 1].recent.end(); itP++) {

						if ((*itP).x + sp > driver[t].end.x) {
							A.x = (*itP).x + sp;
							A.y = (*itP).y - sp;
							driver[t].recent.push_back(A);
						}
					}
					if (A.x != driver[t].end.x) {
						A.x = driver[t].end.x;
						driver[t].recent.push_back(A);
					}

					A.y = driver[t].end.y;
					driver[t].recent.push_back(A);
					//record to dir
					A = driver[t].start;
					driver[t].dir.push_back(A);
					vector <point>::iterator itP;
					for (itP = driver[t].recent.begin();
					itP < driver[t].recent.end(); itP++) {
						driver[t].dir.push_back(*itP);
					}
				}
				k = i + 1;
				a = 1;
			}
		}
	}
	//record to mypath
	for (i = 0; i < driver.size(); i++) {
		a = driver[i].index;
		for (j = 0; j < driver[i].dir.size() - 1; j++) {
			line L;
			L.x1 = driver[i].dir[j].x;
			L.y1 = driver[i].dir[j].y;
			L.x2 = driver[i].dir[j + 1].x;
			L.y2 = driver[i].dir[j + 1].y;
			mypath[a].path.push_back(L);
		}
	}
}

void Walk2(vector < vector <Driver> > &graph, vector<Driver> &driver, vector<Driver> &bump, vector<int> &indexd, vector<int> &indexb, vector <path> &mypath, int width, int space, int wi) {

	int sp = width / 2 + space;//spacing = certain value;
	int a, i, j, k, t, z, z1, ii;
	int position;
	vector<Driver>::iterator itD;

	// determine end;
	a = -1;
	i = 0;
	while (i < driver.size()) {
		z = driver[i].index;
		if ((bump[indexb[z]].pos_b == wi - 2) && (bump[indexb[z]].directlink)) {
			driver[i].end.x = (bump[indexb[z]].x0 + bump[indexb[z]].xf) / 2;
			driver[i].end.y = (bump[indexb[z]].y0 + bump[indexb[z]].yf) / 2;
			a = driver[i].index;// a is the index
			i++;
		}
		else {
			j = i;
			t = 0;
			bool condition = true;
			while (condition) {
				j++;
				t++;
				z1 = driver[j].index;
				if ((bump[indexb[z1]].pos_b == wi - 2) && (bump[indexb[z1]].directlink)) {
					condition = false;
				}
				else if (j == driver.size()) {
					condition = false;
				}
				//record all that goes to the same place as i, there will be t of them and the next is j
			}
			z1 = driver[j].index;
			position = bump[indexb[z1]].pos_a;
			//count back from j
			if (j != driver.size()) {
				for (k = 0; k < t; k++) {
					driver[j - k - 1].end.y = bump[indexb[z1]].yf + (k + 1)* sp;
					vector <point> suitcase;
					suitcase = graph[position - 1][wi - 2].sideleft.interval;
					for (int su = 0; su < suitcase.size(); su++) {
						if (suitcase[su].index == driver[j - k - 1].index) driver[j - k - 1].end.y = suitcase[su].y;
					}
					driver[j - k - 1].end.x = bump[indexb[z1]].xf;
					driver[j - k - 1].end.index = driver[j - k - 1].index;
					graph[position - 1][wi - 2].sideright.interval.push_back(driver[j - k - 1].end);
				}
				i = j;
			}
			else {
				a = driver[i - 1].index;//need to use the previous one;
				position = bump[indexb[a]].pos_a;
				for (k = 0; k < t; k++) {
					driver[i + k].end.y = bump[indexb[a]].y0 - (k + 1)* sp;
					vector <point> suitcase;
					suitcase = graph[position + 1][wi - 2].sideleft.interval;
					for (int su = 0; su < suitcase.size(); su++) {
						if (suitcase[su].index == driver[i + k].index) driver[i + k].end.y = suitcase[su].y;
					}
					driver[i + k].end.x = bump[indexb[a]].xf;
					driver[i + k].end.index = driver[i + k].index;
					graph[position + 1][wi - 2].sideright.interval.push_back(driver[i + k].end);
				}
				i = j;
			}
		}
	}
	//link start to end

	k = 0;
	a = 1;
	for (ii = 0; ii < driver.size(); ii++) {
		if (a == 1) {
			if ((driver[ii].start.y <= driver[ii].end.y) || (ii == driver.size() - 1)) {
				//i-1 往回到k;
				//畫i-1
				//create recent, record to dir;
				if (ii != driver.size() - 1) i = ii - 1;
				else i = ii;
				point A;
				A.y = driver[i].start.y;
				A.x = driver[i].x0 - sp;
				driver[i].recent.push_back(A);
				A.y = driver[i].end.y;
				driver[i].recent.push_back(A);
				A.x = driver[i].end.x;
				driver[i].recent.push_back(A);
				//record to dir		
				A = driver[i].start;
				driver[i].dir.push_back(A);

				vector <point>::iterator itP;
				for (itP = driver[i].recent.begin(); itP < driver[i].recent.end(); itP++) {
					driver[i].dir.push_back(*itP);
				}
				// i-1 to k
				for (t = i - 1; t >= k; t--) {
					A.y = driver[t].start.y;
					A.x = driver[t].x0 - sp;
					driver[t].recent.push_back(A);
					if (driver[t + 1].start.y + sp > driver[t].end.y) {
						A.y = driver[t + 1].start.y + sp;
						driver[t].recent.push_back(A);
					}
					for (itP = driver[t + 1].recent.begin(); itP < driver[t + 1].recent.end(); itP++) {
						if ((*itP).y + sp > driver[t].end.y) {
							A.x = (*itP).x - sp;
							A.y = (*itP).y + sp;
							driver[t].recent.push_back(A);
						}
					}
					if (A.y != driver[t].end.y) {
						A.y = driver[t].end.y;
						driver[t].recent.push_back(A);
					}
					A.x = driver[t].end.x;
					driver[t].recent.push_back(A);
					//record to dir
					A = driver[t].start;
					driver[t].dir.push_back(A);
					vector <point>::iterator itP;
					for (itP = driver[t].recent.begin();
					itP < driver[t].recent.end(); itP++) {
						driver[t].dir.push_back(*itP);

					}
				}
				k = i + 1;
				a = 0;
			}
		}
		else if (a == 0) {
			if ((driver[ii].start.y >= driver[ii].end.y) || (ii == driver.size() - 1)) {
				//k 往後到i;
				//畫k
				//create recent, record to dir;
				if (ii != driver.size() - 1) i = ii - 1;
				else i = ii;
				point A;
				A.y = driver[k].start.y;
				A.x = driver[k].x0 - sp;
				driver[k].recent.push_back(A);
				A.y = driver[k].end.y;
				driver[k].recent.push_back(A);
				A.x = driver[k].end.x;
				driver[k].recent.push_back(A);
				//record to dir
				A = driver[k].start;
				driver[k].dir.push_back(A);
				vector <point>::iterator itP;
				for (itP = driver[k].recent.begin(); itP < driver[k].recent.end(); itP++) {
					driver[k].dir.push_back(*itP);
				}
				//t = k+1 to i
				for (t = k + 1; t <= i; t++) {
					A.y = driver[t].start.y;
					A.x = driver[t].x0 - sp;
					driver[t].recent.push_back(A);
					if (driver[t - 1].start.y - sp < driver[t].end.y) {
						A.y = driver[t - 1].start.y - sp;
						driver[t].recent.push_back(A);
					}
					for (itP = driver[t - 1].recent.begin(); itP < driver[t - 1].recent.end(); itP++) {
						if ((*itP).y - sp < driver[t].end.y) {
							A.x = (*itP).x - sp;
							A.y = (*itP).y - sp;
							driver[t].recent.push_back(A);
						}
					}
					if (A.y != driver[t].end.y) {
						A.y = driver[t].end.y;
						driver[t].recent.push_back(A);
					}
					A.x = driver[t].end.x;
					driver[t].recent.push_back(A);
					//record to dir
					A = driver[t].start;
					driver[t].dir.push_back(A);
					vector <point>::iterator itP;
					for (itP = driver[t].recent.begin();
					itP < driver[t].recent.end(); itP++) {
						driver[t].dir.push_back(*itP);
					}
				}
				k = i + 1;
				a = 1;
			}
		}
	}
	//record to mypath
	for (i = 0; i < driver.size(); i++) {
		a = driver[i].index;
		for (j = 0; j < driver[i].dir.size() - 1; j++) {
			line L;
			L.x1 = driver[i].dir[j].x;
			L.y1 = driver[i].dir[j].y;
			L.x2 = driver[i].dir[j + 1].x;
			L.y2 = driver[i].dir[j + 1].y;
			mypath[a].path.push_back(L);
		}
	}
}

void Walk3(vector < vector <Driver> > &graph, vector<Driver> &driver, vector<Driver> &bump, vector<int> &indexd, vector<int> &indexb, vector <path> &mypath, int width, int space, int le) {

	int sp = width / 2 + space;//spacing = certain value;
	int a, i, j, k, t, z, z1, ii;
	int position;
	vector<Driver>::iterator itD;

	// determine end;
	a = -1;
	i = 0;
	while (i < driver.size()) {
		z = driver[i].index;
		if ((bump[indexb[z]].pos_a == le - 2) && (bump[indexb[z]].directlink)) {
			driver[i].end.x = (bump[indexb[z]].x0 + bump[indexb[z]].xf) / 2;
			driver[i].end.y = (bump[indexb[z]].y0 + bump[indexb[z]].yf) / 2;
			a = driver[i].index;// a is the index
			i++;
		}
		else {
			j = i;
			t = 0;
			bool condition = true;
			while (condition) {
				j++;
				t++;
				z = driver[j].index;
				if ((bump[indexb[z]].pos_a == le - 2) && (bump[indexb[z]].directlink)) {
					condition = false;
				}
				else if (j == driver.size()) {
					condition = false;
				}
				//record same destination as i, #=t next j;
			}
			z1 = driver[j].index;
			position = bump[indexb[z1]].pos_b;
			//back from j
			if (j != driver.size()) {
				for (k = 0; k < t; k++) {
					driver[j - k - 1].end.y = bump[indexb[z1]].y0;
					driver[j - k - 1].end.x = bump[indexb[z1]].xf + (k + 1)* sp;
					vector <point> suitcase;
					suitcase = graph[le - 2][position + 1].sideup.interval;
					for (int su = 0; su < suitcase.size(); su++) {
						if (suitcase[su].index == driver[j - k - 1].index) driver[j - k - 1].end.x = suitcase[su].x;
					}
					driver[j - k - 1].end.index = driver[j - k - 1].index;
					graph[le - 2][position + 1].sidedown.interval.push_back(driver[j - k - 1].end);
				}
				i = j;
			}
			else {
				a = driver[i - 1].index;//need to use the previous one;
				position = bump[indexb[a]].pos_b;
				for (k = 0; k < t; k++) {
					driver[i + k].end.y = bump[indexb[a]].y0;
					driver[i + k].end.x = bump[indexb[a]].x0 - (k + 1)* sp;
					vector <point> suitcase;
					suitcase = graph[le - 2][position - 1].sideup.interval;
					for (int su = 0; su < suitcase.size(); su++) {
						if (suitcase[su].index == driver[i + k].index) driver[i + k].end.x = suitcase[su].x;
					}
					driver[i + k].end.index = driver[i + k].index;
					graph[le - 2][position - 1].sidedown.interval.push_back(driver[i + k].end);
				}
				i = j + 1;
			}
		}
	}

	//link start to end

	k = 0;
	a = 1;
	for (ii = 0; ii < driver.size(); ii++) {
		if (a == 1) {
			if ((driver[ii].start.x <= driver[ii].end.x) || (ii == driver.size() - 1)) {
				//i 往回到k;
				//畫i
				//create recent, record to dir;
				if (ii != driver.size() - 1) i = ii - 1;
				else i = ii;
				point A;
				A.x = driver[i].start.x;
				A.y = driver[i].yf + sp;
				driver[i].recent.push_back(A);
				A.x = driver[i].end.x;
				driver[i].recent.push_back(A);
				A.y = driver[i].end.y;
				driver[i].recent.push_back(A);
				//record to dir		
				A = driver[i].start;
				driver[i].dir.push_back(A);

				vector <point>::iterator itP;
				for (itP = driver[i].recent.begin(); itP < driver[i].recent.end(); itP++) {
					driver[i].dir.push_back(*itP);
				}
				// i-1 to k
				for (t = i - 1; t >= k; t--) {
					A.x = driver[t].start.x;
					A.y = driver[t].yf + sp;
					driver[t].recent.push_back(A);
					if (driver[t + 1].start.x + sp > driver[t].end.x) {
						A.x = driver[t + 1].start.x + sp;
						driver[t].recent.push_back(A);
					}
					for (itP = driver[t + 1].recent.begin(); itP < driver[t + 1].recent.end(); itP++) {
						if ((*itP).x + sp > driver[t].end.x) {
							A.x = (*itP).x + sp;
							A.y = (*itP).y + sp;
							driver[t].recent.push_back(A);
						}
					}
					if (A.x != driver[t].end.x) {
						A.x = driver[t].end.x;
						driver[t].recent.push_back(A);
					}
					A.y = driver[t].end.y;
					driver[t].recent.push_back(A);
					//record to dir
					A = driver[t].start;
					driver[t].dir.push_back(A);
					vector <point>::iterator itP;
					for (itP = driver[t].recent.begin(); itP < driver[t].recent.end(); itP++) {
						driver[t].dir.push_back(*itP);
					}
				}
				k = i + 1;
				a = 0;
			}
		}
		else if (a == 0) {
			if ((driver[ii].start.x >= driver[ii].end.x) || (ii == driver.size() - 1)) {
				//k 往後到i;
				//畫k
				//create recent, record to dir;
				if (ii != driver.size() - 1) i = ii - 1;
				else i = ii;
				point A;
				A.x = driver[k].start.x;
				A.y = driver[k].yf + sp;
				driver[k].recent.push_back(A);
				A.x = driver[k].end.x;
				driver[k].recent.push_back(A);
				A.y = driver[k].end.y;
				driver[k].recent.push_back(A);
				//record to dir
				A = driver[k].start;
				driver[k].dir.push_back(A);
				vector <point>::iterator itP;
				for (itP = driver[k].recent.begin(); itP < driver[k].recent.end(); itP++) {
					driver[k].dir.push_back(*itP);
				}
				//t = k+1 to i
				for (t = k + 1; t <= i; t++) {
					A.x = driver[t].start.x;
					A.y = driver[t].yf + sp;
					driver[t].recent.push_back(A);
					if (driver[t - 1].start.x - sp < driver[t].end.x) {
						A.x = driver[t - 1].start.x - sp;
						driver[t].recent.push_back(A);
					}
					for (itP = driver[t - 1].recent.begin(); itP < driver[t - 1].recent.end(); itP++) {
						if ((*itP).x - sp < driver[t].end.x) {
							A.x = (*itP).x - sp;
							A.y = (*itP).y + sp;
							driver[t].recent.push_back(A);
						}
					}
					if (A.x != driver[t].end.x) {
						A.x = driver[t].end.x;
						driver[t].recent.push_back(A);
					}
					A.y = driver[t].end.y;
					driver[t].recent.push_back(A);
					//record to dir
					A = driver[t].start;
					driver[t].dir.push_back(A);
					vector <point>::iterator itP;
					for (itP = driver[t].recent.begin(); itP < driver[t].recent.end(); itP++) {
						driver[t].dir.push_back(*itP);
					}
				}
				k = i + 1;
				a = 1;
			}
		}
	}

	//record to mypath
	for (i = 0; i < driver.size(); i++) {
		a = driver[i].index;
		for (j = 0; j < driver[i].dir.size() - 1; j++) {
			line L;
			L.x1 = driver[i].dir[j].x;
			L.y1 = driver[i].dir[j].y;
			L.x2 = driver[i].dir[j + 1].x;
			L.y2 = driver[i].dir[j + 1].y;
			mypath[a].path.push_back(L);
		}
	}
}

void Walk4(vector < vector <Driver> > &graph, vector<Driver> &driver, vector<Driver> &bump, vector<int> &indexd, vector<int> &indexb, vector <path> &mypath, int width, int space) {

	int sp = width / 2 + space;//spacing = certain value;
	int a, i, j, k, t, z, z1, ii;
	int position;
	vector<Driver>::iterator itD;

	// determine end;
	a = -1;
	i = 0;
	while (i < driver.size()) {
		z = driver[i].index;
		if ((bump[indexb[z]].pos_b == 1) && (bump[indexb[z]].directlink)) {
			driver[i].end.x = (bump[indexb[z]].x0 + bump[indexb[z]].xf) / 2;
			driver[i].end.y = (bump[indexb[z]].y0 + bump[indexb[z]].yf) / 2;
			a = driver[i].index;// a is the index			
			i++;

		}
		else {
			j = i;
			t = 0;
			bool condition = true;
			while (condition) {
				j++;
				t++;
				z = driver[j].index;
				if ((bump[indexb[z]].pos_b == 1) && (bump[indexb[z]].directlink)) {
					condition = false;
				}
				else if (j == driver.size()) {
					condition = false;
				}
				//record same destination with i, # = t, next is j
			}
			z1 = driver[j].index;
			position = bump[indexb[z1]].pos_a;
			//count back from j
			if (j != driver.size()) {
				z1 = driver[j].index;
				for (k = 0; k < t; k++) {
					driver[j - k - 1].end.y = bump[indexb[z1]].y0 - (k + 1)* sp;
					vector <point> suitcase;
					suitcase = graph[position + 1][1].sideright.interval;
					for (int su = 0; su < suitcase.size(); su++) {
						if (suitcase[su].index == driver[j - k - 1].index) driver[j - k - 1].end.y = suitcase[su].y;
					}
					driver[j - k - 1].end.x = bump[indexb[z1]].x0;
					driver[j - k - 1].end.index = driver[j - k - 1].index;
					graph[position + 1][1].sideleft.interval.push_back(driver[j - k - 1].end);
				}
				i = j;
			}
			else {
				a = driver[i - 1].index;
				position = bump[indexb[a]].pos_a;
				for (k = 0; k < t; k++) {
					driver[i + k].end.y = bump[indexb[a]].yf + (k + 1)* sp;
					vector <point> suitcase;
					suitcase = graph[position - 1][1].sideright.interval;
					for (int su = 0; su < suitcase.size(); su++) {
						if (suitcase[su].index == driver[i + k].index) driver[i + k].end.y = suitcase[su].y;
					}
					driver[i + k].end.x = bump[indexb[a]].x0;
					driver[i + k].end.index = driver[i + k].index;
					graph[position - 1][1].sideleft.interval.push_back(driver[i + k].end);
				}
				i = j;
			}
		}
	}
	//link start to end
	k = 0;
	a = 1;
	for (ii = 0; ii < driver.size(); ii++) {
		if (a == 1) {
			if ((driver[ii].start.y >= driver[ii].end.y) || (ii == driver.size() - 1)) {
				if (ii != driver.size() - 1) i = ii - 1;
				else i = ii;
				//i-1 往回到k;
				//畫i-1
				//create recent, record to dir;
				point A;
				A.y = driver[i].start.y;
				A.x = driver[i].xf + sp;
				driver[i].recent.push_back(A);
				A.y = driver[i].end.y;
				driver[i].recent.push_back(A);
				A.x = driver[i].end.x;
				driver[i].recent.push_back(A);
				//record to dir		
				A = driver[i].start;
				driver[i].dir.push_back(A);
				vector <point>::iterator itP;
				for (itP = driver[i].recent.begin(); itP < driver[i].recent.end(); itP++) {
					driver[i].dir.push_back(*itP);
				}
				// i-1 to k
				for (t = i - 1; t >= k; t--) {
					A.y = driver[t].start.y;
					A.x = driver[t].xf + sp;
					driver[t].recent.push_back(A);
					if (driver[t + 1].start.y - sp < driver[t].end.y) {
						A.y = driver[t + 1].start.y - sp;
						driver[t].recent.push_back(A);
					}
					for (itP = driver[t + 1].recent.begin(); itP < driver[t + 1].recent.end(); itP++) {
						if ((*itP).y - sp < driver[t].end.y) {
							A.x = (*itP).x + sp;
							A.y = (*itP).y - sp;
							driver[t].recent.push_back(A);
						}
					}
					if (A.y != driver[t].end.y) {
						A.y = driver[t].end.y;
						driver[t].recent.push_back(A);
					}
					A.x = driver[t].end.x;
					driver[t].recent.push_back(A);
					//record to dir
					A = driver[t].start;
					driver[t].dir.push_back(A);
					vector <point>::iterator itP;
					for (itP = driver[t].recent.begin(); itP < driver[t].recent.end(); itP++) {
						driver[t].dir.push_back(*itP);

					}
				}
				k = i + 1;
				a = 0;
			}
		}
		else if (a == 0) {
			if ((driver[ii].start.y <= driver[ii].end.y) || (ii == driver.size() - 1)) {
				//k 往後到i;
				//畫k
				//create recent, record to dir;				
				if (ii != driver.size() - 1) i = ii - 1;
				else i = ii;
				point A;
				A.y = driver[k].start.y;
				A.x = driver[k].xf + sp;
				driver[k].recent.push_back(A);
				A.y = driver[k].end.y;
				driver[k].recent.push_back(A);
				A.x = driver[k].end.x;
				driver[k].recent.push_back(A);
				//record to dir
				A = driver[k].start;
				driver[k].dir.push_back(A);
				vector <point>::iterator itP;
				for (itP = driver[k].recent.begin(); itP < driver[k].recent.end(); itP++) {
					driver[k].dir.push_back(*itP);
				}
				//t = k+1 to i
				for (t = k + 1; t <= i; t++) {
					A.y = driver[t].start.y;
					A.x = driver[t].xf + sp;
					driver[t].recent.push_back(A);
					if (driver[t - 1].start.y + sp > driver[t].end.y) {
						A.y = driver[t - 1].start.y + sp;
						driver[t].recent.push_back(A);
					}
					for (itP = driver[t - 1].recent.begin(); itP < driver[t - 1].recent.end(); itP++) {
						if ((*itP).y + sp > driver[t].end.y) {
							A.x = (*itP).x + sp;
							A.y = (*itP).y + sp;
							driver[t].recent.push_back(A);
						}
					}
					if (A.y != driver[t].end.y) {
						A.y = driver[t].end.y;
						driver[t].recent.push_back(A);
					}
					A.x = driver[t].end.x;
					driver[t].recent.push_back(A);
					//record to dir
					A = driver[t].start;
					driver[t].dir.push_back(A);
					vector <point>::iterator itP;
					for (itP = driver[t].recent.begin(); itP < driver[t].recent.end(); itP++) {
						driver[t].dir.push_back(*itP);
					}
				}
				k = i + 1;
				a = 1;
			}
		}
	}

	//record to mypath
	for (i = 0; i < driver.size(); i++) {
		a = driver[i].index;
		for (j = 0; j < driver[i].dir.size() - 1; j++) {
			line L;
			L.x1 = driver[i].dir[j].x;
			L.y1 = driver[i].dir[j].y;
			L.x2 = driver[i].dir[j + 1].x;
			L.y2 = driver[i].dir[j + 1].y;
			mypath[a].path.push_back(L);
		}
	}
}

void Sort_line(vector <detail> &line, string order);

void fill(vector <vector <Driver> > &graph,int i,int j){

  	point node;

  	for(int k=0;k<graph[i-1][j].up.size();k++){     //above upward
    	node=graph[i-1][j].up[k];
    	graph[i][j].sideup.outside.push_back(node);
  	}
  	for(int k=0;k<graph[i-1][j].down.size();k++){   //above downward
    	node=graph[i-1][j].down[k];
    	graph[i][j].sideup.inside.push_back(node);
  	}	
  	for(int k=0;k<graph[i][j+1].right.size();k++){  //right to right
    	node=graph[i][j+1].right[k];
    	graph[i][j].sideright.outside.push_back(node);
  	}
  	for(int k=0;k<graph[i][j+1].left.size();k++){    //right to left
    	node=graph[i][j+1].left[k];
    	graph[i][j].sideright.inside.push_back(node);
  	}		
	for(int k=0;k<graph[i+1][j].down.size();k++){  //down  upward
    	node=graph[i+1][j].down[k];
    	graph[i][j].sidedown.outside.push_back(node);
  	}
  	for(int k=0;k<graph[i+1][j].up.size();k++){    //down downward
    	node=graph[i+1][j].up[k];
    	graph[i][j].sidedown.inside.push_back(node);
  	}
	for(int k=0;k<graph[i][j-1].left.size();k++){  // left to left
    	node=graph[i][j-1].left[k];
    	graph[i][j].sideleft.outside.push_back(node);
  	}
  	for(int k=0;k<graph[i][j-1].right.size();k++){    //left to right
    	node=graph[i][j-1].right[k];
    	graph[i][j].sideleft.inside.push_back(node);
  	}                             
}

void detailed(vector <vector <Driver> > &graph, vector < path > &mypath ,int i,int j,int width,int space){
	
	vector <detail> LU;
	vector <detail> RU;
	vector <detail> RD;
	vector <detail> LD;
	vector <detail> LR;
	vector <detail> UD;
	vector <detail> UL;
	vector <detail> UR;
	vector <detail> DR;
	vector <detail> DL;
	vector <detail> RL;
	vector <detail> DU;
	vector <detail> Single;

	vector <detail> out;
	vector <detail> in; 

	//assign all node into vector out,in
	for(int k=0;k<graph[i][j].sideup.outside.size();k++){  
		detail tmp;
		tmp.sec=1;
		tmp.index=graph[i][j].sideup.outside[k].index;
		tmp.virtual_index=graph[i][j].sideup.outside[k].virtual_index;
		out.push_back(tmp);
	}

	for(int k=0;k<graph[i][j].sideright.outside.size();k++){
		detail tmp;
		tmp.sec=2;
		tmp.index=graph[i][j].sideright.outside[k].index;
		tmp.virtual_index=graph[i][j].sideright.outside[k].virtual_index;
		out.push_back(tmp);
	}

	for(int k=0;k<graph[i][j].sidedown.outside.size();k++){
		detail tmp;
		tmp.sec=3;
		tmp.index=graph[i][j].sidedown.outside[k].index;
		tmp.virtual_index=graph[i][j].sidedown.outside[k].virtual_index;
		out.push_back(tmp);
	}

	for(int k=0;k<graph[i][j].sideleft.outside.size();k++){
		detail tmp;
		tmp.sec=4;
		tmp.index=graph[i][j].sideleft.outside[k].index;
		tmp.virtual_index=graph[i][j].sideleft.outside[k].virtual_index;
		out.push_back(tmp);
	}

	for(int k=0;k<graph[i][j].sideup.inside.size();k++){
		detail tmp;
		tmp.sec=1;
		tmp.index=graph[i][j].sideup.inside[k].index;
		tmp.virtual_index=graph[i][j].sideup.inside[k].virtual_index;
		in.push_back(tmp);
	}

	for(int k=0;k<graph[i][j].sideright.inside.size();k++){
		detail tmp;
		tmp.sec=2;
		tmp.index=graph[i][j].sideright.inside[k].index;
		tmp.virtual_index=graph[i][j].sideright.inside[k].virtual_index;
		in.push_back(tmp);
	}

	for(int k=0;k<graph[i][j].sidedown.inside.size();k++){
		detail tmp;
		tmp.sec=3;
		tmp.index=graph[i][j].sidedown.inside[k].index;
		tmp.virtual_index=graph[i][j].sidedown.inside[k].virtual_index;
		in.push_back(tmp);
	}

	for(int k=0;k<graph[i][j].sideleft.inside.size();k++){
		detail tmp;
		tmp.sec=4;
		tmp.index=graph[i][j].sideleft.inside[k].index;
		tmp.virtual_index=graph[i][j].sideleft.inside[k].virtual_index;
		in.push_back(tmp);
	}
	//END assign all node into vector out,in

	vector <int> pillar(4);
	pillar[0]=graph[i-1][j-1].index;
	pillar[1]=graph[i-1][j+1].index;
	pillar[2]=graph[i+1][j+1].index;
	pillar[3]=graph[i+1][j-1].index;
	int x0=graph[i][j].x0;
	int y0=graph[i][j].y0;
	int xf=graph[i][j].xf;
	int yf=graph[i][j].yf;
	int interval=space+width;
	int sec1_lx=x0+space+width/2;
	int sec1_rx=xf-space-width/2;
	int sec2_uy=yf-space-width/2;
	int sec2_dy=y0+space+width/2;
	int sec3_lx=x0+space+width/2;
	int sec3_rx=xf-space-width/2;
	int sec4_uy=yf-space-width/2;
	int sec4_dy=y0+space+width/2;
	
	//find all node that only appear once(only in or only out)
		for(int k=0;k<in.size();k++){
			detail ind=in[k];
			int count=0;
			for(int l=0;l<out.size();l++){
				if(out[l].virtual_index==ind.virtual_index)
					count++;
			}
			if (count==0){
				Single.push_back(ind);
			}
		}
		
		for(int k=0;k < out.size() ; k++ ){
			detail ind=out[k];
			int count=0;
			for(int l=0;l<in.size();l++){
				if(in[l].virtual_index==ind.virtual_index)
					count++;
			}
			if (count==0)
				Single.push_back(ind);
		}
	//END find all node that only appear once(only in or only out)
	
	//Start decide the node which only in or out
	point tmp_point;
	for(int k=0;k<Single.size();k++){      
		if(Single[k].index ==	pillar[0]){
			if(Single[k].sec==1){ 
				Single[k].x=sec1_lx-space;
				Single[k].y=yf;
				tmp_point.x=sec1_lx-space;
				tmp_point.y=yf;
				tmp_point.index=pillar[0];
				tmp_point.virtual_index=Single[k].virtual_index;
				graph[i-1][j].sidedown.interval.push_back(tmp_point);
				sec1_lx+=interval;
			}
			if(Single[k].sec==4){
				Single[k].x=x0;
				Single[k].y=sec4_uy+space;
				tmp_point.x=x0;
				tmp_point.y=sec4_uy+space;
				tmp_point.index=pillar[0];
				tmp_point.virtual_index=Single[k].virtual_index;
				graph[i][j-1].sideright.interval.push_back(tmp_point);
				sec4_uy-=interval;
			}
		}
		if(Single[k].index ==	pillar[1]){
			if(Single[k].sec==1){ 
				Single[k].x=sec1_rx+space;
				Single[k].y=yf;
				tmp_point.x=sec1_rx+space;
				tmp_point.y=yf;
				tmp_point.index=pillar[1];
				tmp_point.virtual_index=Single[k].virtual_index;
				graph[i-1][j].sidedown.interval.push_back(tmp_point);
				sec1_rx-=interval;
			}
			if(Single[k].sec==2){
				Single[k].x=xf;
				Single[k].y=sec2_uy+space;
				tmp_point.x=xf;
				tmp_point.y=sec2_uy+space;
				tmp_point.index=pillar[1];
				tmp_point.virtual_index=Single[k].virtual_index;
				graph[i][j+1].sideleft.interval.push_back(tmp_point);
				sec2_uy-=interval;
			}
		}
		if(Single[k].index ==	pillar[2]){
			if(Single[k].sec==2){ 
				Single[k].x=xf;
				Single[k].y=sec2_dy-space;
				tmp_point.x=xf;
				tmp_point.y=sec2_dy-space;
				tmp_point.index=pillar[2];
				tmp_point.virtual_index=Single[k].virtual_index;
				graph[i][j+1].sideleft.interval.push_back(tmp_point);
				sec2_dy+=interval;
			}
			if(Single[k].sec==3){
				Single[k].x=sec3_rx+space;
				Single[k].y=y0;
				tmp_point.x=sec3_rx+space;
				tmp_point.y=y0;
				tmp_point.index=pillar[2];
				tmp_point.virtual_index=Single[k].virtual_index;
				graph[i+1][j].sideup.interval.push_back(tmp_point);
				sec3_rx-=interval;
			}
		}
		if(Single[k].index ==	pillar[3]){
			if(Single[k].sec==3){ 
				Single[k].x=sec3_lx-space;
				Single[k].y=y0;
				tmp_point.x=sec3_lx-space;
				tmp_point.y=y0;
				tmp_point.index=pillar[3];
				tmp_point.virtual_index=Single[k].virtual_index;
				graph[i+1][j].sideup.interval.push_back(tmp_point);
				sec3_lx+=interval;
			}
			if(Single[k].sec==4){
				Single[k].x=x0;
				Single[k].y=sec4_dy-space;
				tmp_point.x=x0;
				tmp_point.y=sec4_dy-space;
				tmp_point.index=pillar[3];
				tmp_point.virtual_index=Single[k].virtual_index;
				graph[i][j-1].sideright.interval.push_back(tmp_point);
				sec4_dy+=interval;
			}
		}
	}
	//END decide the node which only in or out

	//assign all path to 12 type
	for(int k=0;k<in.size();k++){
		for(int l=0;l<out.size();l++){
			if(in[k].virtual_index==out[l].virtual_index){
				if(in[k].sec==4 && out[l].sec==1)
					LU.push_back(in[k]);
				if(in[k].sec==1 && out[l].sec==2)
					UR.push_back(in[k]);
				if(in[k].sec==2 && out[l].sec==3)
					RD.push_back(in[k]);
				if(in[k].sec==3 && out[l].sec==4)
					DL.push_back(in[k]);
				if(in[k].sec==1 && out[l].sec==3)
					UD.push_back(in[k]);
				if(in[k].sec==4 && out[l].sec==2)
					LR.push_back(in[k]);	
				if(in[k].sec==1 && out[l].sec==4)
					UL.push_back(in[k]);
				if(in[k].sec==2 && out[l].sec==1)
					RU.push_back(in[k]);
				if(in[k].sec==3 && out[l].sec==2)
					DR.push_back(in[k]);
				if(in[k].sec==4 && out[l].sec==3)
					LD.push_back(in[k]);
				if(in[k].sec==3 && out[l].sec==1)
					DU.push_back(in[k]);
				if(in[k].sec==2 && out[l].sec==4)
					RL.push_back(in[k]);	
			}


		}
	}
	//END assign all path to 12 type

	//Sort according to its direction
	Sort_line(LU,"O");
	Sort_line(UR,"O");
	Sort_line(RD,"O");
	Sort_line(DL,"O");
	Sort_line(UL,"U");
	Sort_line(RU,"U");
	Sort_line(DR,"U");
	Sort_line(LD,"U");
	Sort_line(UD,"U");
	Sort_line(DU,"U");
	Sort_line(LR,"U");
	Sort_line(RL,"U");
	//End Sort according to its direction

	//assign  coordinate to the path
	line tmp;
	for(int k=0;k<LU.size();k++){
		tmp.x1=sec1_lx;
		tmp.x2=sec1_lx;
		tmp.y1=yf;
		tmp.y2=sec4_uy;
		int ind=LU[k].index;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=sec1_lx;
		tmp_point.y=yf;
		tmp_point.index=ind;
		graph[i-1][j].sidedown.interval.push_back(tmp_point);

		tmp.x1=x0;
		tmp.x2=sec1_lx;
		tmp.y1=sec4_uy;
		tmp.y2=sec4_uy;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=x0;
		tmp_point.y=sec4_uy;
		tmp_point.index=ind;
		graph[i][j-1].sideright.interval.push_back(tmp_point);

		sec1_lx+=interval;
		sec4_uy-=interval;	
	}
	
	for(int k=0;k<UL.size();k++){
		tmp.x1=sec1_lx;
		tmp.x2=sec1_lx;
		tmp.y1=yf;
		tmp.y2=sec4_uy;
		int ind=UL[k].index;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=sec1_lx;
		tmp_point.y=yf;
		tmp_point.index=ind;
		graph[i-1][j].sidedown.interval.push_back(tmp_point);

		tmp.x1=x0;
		tmp.x2=sec1_lx;
		tmp.y1=sec4_uy;
		tmp.y2=sec4_uy;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=x0;
		tmp_point.y=sec4_uy;
		tmp_point.index=ind;
		graph[i][j-1].sideright.interval.push_back(tmp_point);

		sec1_lx+=interval;
		sec4_uy-=interval;		
	}
	
	for(int k=0;k<UR.size();k++){
		tmp.x1=sec1_rx;
		tmp.x2=sec1_rx;
		tmp.y1=yf;
		tmp.y2=sec2_uy;
		int ind=UR[k].index;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=sec1_rx;
		tmp_point.y=yf;
		tmp_point.index=ind;
		graph[i-1][j].sidedown.interval.push_back(tmp_point);

		tmp.x1=xf;
		tmp.x2=sec1_rx;
		tmp.y1=sec2_uy;
		tmp.y2=sec2_uy;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=xf;
		tmp_point.y=sec2_uy;
		tmp_point.index=ind;
		graph[i][j+1].sideleft.interval.push_back(tmp_point);

		sec1_rx-=interval;
		sec2_uy-=interval;
	}
	
	for(int k=0;k<RU.size();k++){
		tmp.x1=sec1_rx;
		tmp.x2=sec1_rx;
		tmp.y1=yf;
		tmp.y2=sec2_uy;
		int ind=RU[k].index;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=sec1_rx;
		tmp_point.y=yf;
		tmp_point.index=ind;
		graph[i-1][j].sidedown.interval.push_back(tmp_point);

		tmp.x1=xf;
		tmp.x2=sec1_rx;
		tmp.y1=sec2_uy;
		tmp.y2=sec2_uy;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=xf;
		tmp_point.y=sec2_uy;
		tmp_point.index=ind;
		graph[i][j+1].sideleft.interval.push_back(tmp_point);

		sec1_rx-=interval;
		sec2_uy-=interval;
	}
	
	for(int k=0;k<RD.size();k++){
		tmp.x1=xf;
		tmp.x2=sec3_rx;
		tmp.y1=sec2_dy;
		tmp.y2=sec2_dy;
		int ind=RD[k].index;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=xf;
		tmp_point.y=sec2_dy;
		tmp_point.index=ind;
		graph[i][j+1].sideleft.interval.push_back(tmp_point);

		tmp.x1=sec3_rx;
		tmp.x2=sec3_rx;
		tmp.y1=y0;
		tmp.y2=sec2_dy;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=sec3_rx;
		tmp_point.y=y0;
		tmp_point.index=ind;
		graph[i+1][j].sideup.interval.push_back(tmp_point);

		sec3_rx-=interval;
		sec2_dy+=interval;	
	}
	
	for(int k=0;k<DR.size();k++){
		tmp.x1=xf;
		tmp.x2=sec3_rx;
		tmp.y1=sec2_dy;
		tmp.y2=sec2_dy;
		int ind=DR[k].index;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=xf;
		tmp_point.y=sec2_dy;
		tmp_point.index=ind;
		graph[i][j+1].sideleft.interval.push_back(tmp_point);

		tmp.x1=sec3_rx;
		tmp.x2=sec3_rx;
		tmp.y1=y0;
		tmp.y2=sec2_dy;
		mypath[ind].path.push_back(tmp);		
		tmp_point.x=sec3_rx;
		tmp_point.y=y0;
		tmp_point.index=ind;
		graph[i+1][j].sideup.interval.push_back(tmp_point);

		sec3_rx-=interval;
		sec2_dy+=interval;	
	}

	for(int k=0;k<DL.size();k++){
		tmp.x1=sec3_lx;
		tmp.x2=sec3_lx;
		tmp.y1=y0;
		tmp.y2=sec4_dy;
		int ind=DL[k].index;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=sec3_lx;
		tmp_point.y=y0;
		tmp_point.index=ind;
		graph[i+1][j].sideup.interval.push_back(tmp_point);

		tmp.x1=x0;
		tmp.x2=sec3_lx;
		tmp.y1=sec4_dy;
		tmp.y2=sec4_dy;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=x0;
		tmp_point.y=sec4_dy;
		tmp_point.index=ind;
		graph[i][j-1].sideright.interval.push_back(tmp_point);

		sec3_lx+=interval;
		sec4_dy+=interval;	
	}

	for(int k=0;k<LD.size();k++){
		tmp.x1=sec3_lx;
		tmp.x2=sec3_lx;
		tmp.y1=y0;
		tmp.y2=sec4_dy;
		int ind=LD[k].index;
		mypath[ind].path.push_back(tmp);
		tmp_point.x=sec3_lx;
		tmp_point.y=y0;
		tmp_point.index=ind;
		graph[i+1][j].sideup.interval.push_back(tmp_point);

		tmp.x1=x0;
		tmp.x2=sec3_lx;
		tmp.y1=sec4_dy;
		tmp.y2=sec4_dy;
		mypath[ind].path.push_back(tmp);						
		tmp_point.x=x0;
		tmp_point.y=sec4_dy;
		tmp_point.index=ind;
		graph[i][j-1].sideright.interval.push_back(tmp_point);

		sec3_lx+=interval;
		sec4_dy+=interval;	
	}

	for(int k=0;k < DU.size();k++){
		tmp_point.x=sec3_rx;
		tmp_point.y=y0;
		tmp_point.index=DU[k].index;
		graph[i+1][j].sideup.interval.push_back(tmp_point);
		tmp_point.x=sec1_rx;
		tmp_point.y=yf;
		tmp_point.index=DU[k].index;
		graph[i-1][j].sidedown.interval.push_back(tmp_point);
		
		if(sec3_rx==sec1_rx){
			tmp.x1=sec3_rx;
			tmp.x2=sec3_rx;
			tmp.y1=y0;
			tmp.y2=yf;
			mypath[DU[k].index].path.push_back(tmp);


			sec3_rx-=interval;
			sec1_rx-=interval;
		}
		else{
			if(sec1_rx > sec3_rx){
				tmp.x1=sec3_rx;
				tmp.x2=sec3_rx;
				tmp.y1=y0;
				tmp.y2=sec2_dy;
				mypath[DU[k].index].path.push_back(tmp);
				tmp.x1=sec3_rx;
				tmp.x2=sec1_rx;
				tmp.y1=sec2_dy;
				tmp.y2=sec2_dy;
				mypath[DU[k].index].path.push_back(tmp);
				tmp.x1=sec1_rx;
				tmp.x2=sec1_rx;
				tmp.y1=sec2_dy;
				tmp.y2=yf;
				mypath[DU[k].index].path.push_back(tmp);
				sec1_rx-=interval;
				sec2_dy+=interval;
				sec3_rx-=interval;
			}
			else{
				tmp.x1=sec3_rx;
				tmp.x2=sec3_rx;
				tmp.y1=y0;
				tmp.y2=sec2_uy;
				mypath[DU[k].index].path.push_back(tmp);
				tmp.x1=sec3_rx;
				tmp.x2=sec1_rx;
				tmp.y1=sec2_uy;
				tmp.y2=sec2_uy;
				mypath[DU[k].index].path.push_back(tmp);
				tmp.x1=sec1_rx;
				tmp.x2=sec1_rx;
				tmp.y1=sec2_uy;
				tmp.y2=yf;
				mypath[DU[k].index].path.push_back(tmp);				
				sec1_rx-=interval;
				sec3_rx-=interval;
				sec2_uy-=interval;
			}
		}
	}
	
	for(int k=0;k < UD.size();k++){
		tmp_point.x=sec3_lx;
		tmp_point.y=y0;
		tmp_point.index=UD[k].index;
		graph[i+1][j].sideup.interval.push_back(tmp_point);
		tmp_point.x=sec1_lx;
		tmp_point.y=yf;
		tmp_point.index=UD[k].index;
		graph[i-1][j].sidedown.interval.push_back(tmp_point);
	
		if(sec3_lx==sec1_lx){
			tmp.x1=sec3_lx;
			tmp.x2=sec3_lx;
			tmp.y1=y0;
			tmp.y2=yf;
			mypath[UD[k].index].path.push_back(tmp);
			sec3_lx+=interval;
			sec1_lx+=interval;
		}
		else{
			if(sec1_lx < sec3_lx){
				tmp.x1=sec3_lx;
				tmp.x2=sec3_lx;
				tmp.y1=y0;
				tmp.y2=sec4_dy;
				mypath[UD[k].index].path.push_back(tmp);
				tmp.x1=sec3_lx;
				tmp.x2=sec1_lx;
				tmp.y1=sec4_dy;
				tmp.y2=sec4_dy;
				mypath[UD[k].index].path.push_back(tmp);
				tmp.x1=sec1_lx;
				tmp.x2=sec1_lx;
				tmp.y1=sec4_dy;
				tmp.y2=yf;
				mypath[UD[k].index].path.push_back(tmp);
				sec1_lx+=interval;
				sec4_dy+=interval;
				sec3_lx+=interval;
			}
			else{
				tmp.x1=sec3_lx;
				tmp.x2=sec3_lx;
				tmp.y1=y0;
				tmp.y2=sec4_uy;
				mypath[UD[k].index].path.push_back(tmp);
				tmp.x1=sec3_lx;
				tmp.x2=sec1_lx;
				tmp.y1=sec4_uy;
				tmp.y2=sec4_uy;
				mypath[UD[k].index].path.push_back(tmp);
				tmp.x1=sec1_lx;
				tmp.x2=sec1_lx;
				tmp.y1=sec4_uy;
				tmp.y2=yf;
				mypath[UD[k].index].path.push_back(tmp);				
				sec1_lx+=interval;
				sec4_uy-=interval;
				sec3_lx+=interval;
			}
		}
	}
	
	for(int k=0;k < LR.size();k++){
		tmp_point.x=x0;
		tmp_point.y=sec4_dy;
		tmp_point.index=LR[k].index;
		graph[i][j-1].sideright.interval.push_back(tmp_point);
		tmp_point.x=xf;
		tmp_point.y=sec2_dy;
		tmp_point.index=LR[k].index;
		graph[i][j+1].sideleft.interval.push_back(tmp_point);
		
		if(sec4_dy==sec2_dy){
			tmp.x1=x0;
			tmp.x2=xf;
			tmp.y1=sec4_dy;
			tmp.y2=sec4_dy;
			mypath[LR[k].index].path.push_back(tmp);
			sec4_dy+=interval;
			sec2_dy+=interval;
		}
		else{
			if(sec4_dy > sec2_dy){
				tmp.x1=x0;
				tmp.x2=sec3_lx;
				tmp.y1=sec4_dy;
				tmp.y2=sec4_dy;
				mypath[LR[k].index].path.push_back(tmp);
				tmp.x1=sec3_lx;
				tmp.x2=sec3_lx;
				tmp.y1=sec4_dy;
				tmp.y2=sec2_dy;
				mypath[LR[k].index].path.push_back(tmp);
				tmp.x1=sec3_lx;
				tmp.x2=xf;
				tmp.y1=sec2_dy;
				tmp.y2=sec2_dy;
				mypath[LR[k].index].path.push_back(tmp);
				sec4_dy+=interval;
				sec2_dy+=interval;
				sec3_lx+=interval;
			}
			else{
				tmp.x1=x0;
				tmp.x2=sec3_rx;
				tmp.y1=sec4_dy;
				tmp.y2=sec4_dy;
				mypath[LR[k].index].path.push_back(tmp);
				tmp.x1=sec3_rx;
				tmp.x2=sec3_rx;
				tmp.y1=sec4_dy;
				tmp.y2=sec2_dy;
				mypath[LR[k].index].path.push_back(tmp);
				tmp.x1=sec3_rx;
				tmp.x2=xf;
				tmp.y1=sec2_dy;
				tmp.y2=sec2_dy;
				mypath[LR[k].index].path.push_back(tmp);				
				sec4_dy+=interval;
				sec2_dy+=interval;
				sec3_rx-=interval;
			}
		}
	}

	for(int k=0;k < RL.size();k++){
		tmp_point.x=x0;
		tmp_point.y=sec4_uy;
		tmp_point.index=RL[k].index;
		graph[i][j-1].sideright.interval.push_back(tmp_point);
		tmp_point.x=xf;
		tmp_point.y=sec2_uy;
		tmp_point.index=RL[k].index;
		graph[i][j+1].sideleft.interval.push_back(tmp_point);
		
		if(sec4_uy==sec2_uy){
			tmp.x1=xf;
			tmp.x2=x0;
			tmp.y1=sec4_uy;
			tmp.y2=sec4_uy;
			mypath[RL[k].index].path.push_back(tmp);
			sec4_uy-=interval;
			sec2_uy-=interval;
		}
		else{
			if(sec2_uy < sec4_uy){
				tmp.x1=xf;
				tmp.x2=sec1_rx;
				tmp.y1=sec2_uy;
				tmp.y2=sec2_uy;
				mypath[RL[k].index].path.push_back(tmp);
				tmp.x1=sec1_rx;
				tmp.x2=sec1_rx;
				tmp.y1=sec2_uy;
				tmp.y2=sec4_uy;
				mypath[RL[k].index].path.push_back(tmp);
				tmp.x1=sec1_rx;
				tmp.x2=x0;
				tmp.y1=sec4_uy;
				tmp.y2=sec4_uy;
				mypath[RL[k].index].path.push_back(tmp);
				sec4_uy-=interval;
				sec2_uy-=interval;
				sec1_rx-=interval;
			}
			else{
				tmp.x1=xf;
				tmp.x2=sec1_lx;
				tmp.y1=sec2_uy;
				tmp.y2=sec2_uy;
				mypath[RL[k].index].path.push_back(tmp);
				tmp.x1=sec1_lx;
				tmp.x2=sec1_lx;
				tmp.y1=sec2_uy;
				tmp.y2=sec4_uy;
				mypath[RL[k].index].path.push_back(tmp);
				tmp.x1=sec1_lx;
				tmp.x2=x0;
				tmp.y1=sec4_uy;
				tmp.y2=sec4_uy;
				mypath[RL[k].index].path.push_back(tmp);				
				sec4_uy-=interval;
				sec2_uy-=interval;
				sec1_lx+=interval;
			}
		}
	}
	//End assign  coordinate to the path
}

void Sort_line(vector <detail> &line, string order){
	if (order=="O"){
		for(int i=0;i<line.size();i++){    
			detail tmp=line[i];   					
			int j=i-1;
			while ( j>-1 && (tmp.virtual_index < line[j].virtual_index)){
				line[j+1]=line[j];
				line[j]=tmp;
				j--;
			}
		}
	}
	else{
		for(int i=0;i<line.size();i++){    
			detail tmp=line[i];   					
			int j=i-1;
			while ( j>-1 && (tmp.virtual_index > line[j].virtual_index)){
				line[j+1]=line[j];
				line[j]=tmp;
				j--;
			}
		}
	}		
}

void Sort_y(vector<Driver> &bump, string order){
	if (order=="O"){
		for(int i=0;i<bump.size();i++){    
			Driver tmp=bump[i];   					
			int j=i-1;
			while ( j>-1 && (tmp.yf < bump[j].yf)){
				bump[j+1]=bump[j];
				bump[j]=tmp;
				j--;
			}
		}
	}
	else{
		for(int i=0;i<bump.size();i++){    
			Driver tmp=bump[i];   					
			int j=i-1;
			while ( j>-1 && (tmp.yf > bump[j].yf)){
				bump[j+1]=bump[j];
				bump[j]=tmp;
				j--;
			}
		}
	}
}

void Sort_x(vector<Driver> &bump, string order){
	if (order=="O"){
		for(int i=0;i<bump.size();i++){    
			Driver tmp=bump[i];   					
			int j=i-1;
			while ( j>-1 && (tmp.xf < bump[j].xf)){
				bump[j+1]=bump[j];
				bump[j]=tmp;
				j--;
			}
		}
	}
	else{
		for(int i=0;i<bump.size();i++){    
			Driver tmp=bump[i];   					
			int j=i-1;
			while ( j>-1 && (tmp.xf > bump[j].xf)){
				bump[j+1]=bump[j];
				bump[j]=tmp;
				j--;
			}
		}
	}
}

int SetRing(vector<Driver> &bump, string driver_section, vector<Driver> &cross_virtual_bump){

	int ring_index=1;
	int single_ring_num=1;
	bump[0].ring=ring_index;
	bump[0].index_in_ring=single_ring_num;
	if(driver_section=="L" || driver_section=="R"){
		for(int i=1;i<bump.size();i++){
			if(!(bump[i].xf == bump[i-1].xf)){
				bump[i].ring=(++ring_index);
				single_ring_num=1;
				bump[i].index_in_ring=single_ring_num;
			}
			else {
				bump[i].ring=ring_index;
				bump[i].index_in_ring=(++single_ring_num);	
			}
		}
	}
	else{
		for(int i=1;i<bump.size();i++){
			if(!(bump[i].yf == bump[i-1].yf)){
				bump[i].ring=(++ring_index);
				single_ring_num=1;
				bump[i].index_in_ring=single_ring_num;
			}
			else{
				bump[i].ring=ring_index;
				bump[i].index_in_ring=(++single_ring_num);	
			}
		}
	}
	ring_index++;
	single_ring_num=1;
	int m = bump.size()-1;
	while (bump[m].ring == ring_index-1){
		m--;
	}
	if (driver_section=="U"){
		int n = bump[m+1].pos_a + 2;
		int o = bump[m+1].pos_b;
		int size=bump.size();
		for(int i=0;i< size;i++){
			if(bump[i].sec[0]!=bump[i].sec[1]){
				int j=i;
				bump[i].real=false;
				while(bump[i].ring == bump[j].ring && j<=size){
					if(bump[i].virtual_index >bump[j].index){
						bump[i].virtual_index=bump[j].index-0.75;
					}
					j++;
				}
                Driver tmp;
                tmp.index=bump[i].index;
                tmp.sec.push_back(bump[i].sec[0]);
                tmp.sec.push_back(bump[i].sec[1]);
                tmp.virtual_index=bump[i].virtual_index+0.25;//8.5 78.25
                tmp.pos_a=n;
                tmp.pos_b=o;
                tmp.ring=ring_index;
                tmp.index_in_ring=single_ring_num;
                o+=2;
                bump.push_back(tmp);
                cross_virtual_bump.push_back(tmp);
                single_ring_num++;
			}
		}
	}
	else if (driver_section=="R"){
		int n = bump[m+1].pos_a;
		int o = bump[m+1].pos_b-2;
		int size=bump.size();
		for(int i=0;i< size;i++){
			if(bump[i].sec[0]!=bump[i].sec[1]){
				int j=i;
				bump[i].real=false;
				while(bump[i].ring == bump[j].ring && j<=size){
					if(bump[i].virtual_index >bump[j].index){
						bump[i].virtual_index=bump[j].index-0.75;
					}
					j++;
				}
			Driver tmp;
            tmp.sec.push_back(bump[i].sec[0]);
            tmp.sec.push_back(bump[i].sec[1]);
            tmp.index=bump[i].index;
			tmp.virtual_index=bump[i].virtual_index+0.25;
			tmp.pos_a=n;
			tmp.pos_b=o;
			tmp.ring=ring_index;
			tmp.index_in_ring=single_ring_num;
			n+=2;
			bump.push_back(tmp);
			cross_virtual_bump.push_back(tmp);
			single_ring_num++;
			}
		}
	}
	else if (driver_section=="D"){
		int n = bump[m+1].pos_a-2;
		int o = bump[m+1].pos_b;
		int size=bump.size();
		for(int i=0;i< size;i++){
			if(bump[i].sec[0]!=bump[i].sec[1]){
				int j=i;
				bump[i].real=false;
				while(bump[i].ring == bump[j].ring && j<=size){
					if(bump[i].virtual_index >bump[j].index){
						bump[i].virtual_index=bump[j].index-0.75;
					}
					j++;
				}
			Driver tmp;
            tmp.sec.push_back(bump[i].sec[0]);
            tmp.sec.push_back(bump[i].sec[1]);
            tmp.index=bump[i].index;
			tmp.virtual_index=bump[i].virtual_index+0.25;
			tmp.pos_a=n;
			tmp.pos_b=o;
			tmp.ring=ring_index;
			tmp.index_in_ring=single_ring_num;
			o+=2;
			bump.push_back(tmp);
			cross_virtual_bump.push_back(tmp);
			single_ring_num++;
			}
		}
	}
	else if (driver_section=="L"){
		int n = bump[m+1].pos_a;
		int o = bump[m+1].pos_b+2;
		int size=bump.size();
		for(int i=0;i< size;i++){
			if(bump[i].sec[0]!=bump[i].sec[1]){
				int j=i;
				bump[i].real=false;
				while(bump[i].ring == bump[j].ring && j<=size){
					if(bump[i].virtual_index >bump[j].index){
						bump[i].virtual_index=bump[j].index-0.75;
					}
					j++;
				}
			Driver tmp;
            tmp.sec.push_back(bump[i].sec[0]);
            tmp.sec.push_back(bump[i].sec[1]);
            tmp.index=bump[i].index;
			tmp.virtual_index=bump[i].virtual_index+0.25;
			tmp.pos_a=n;
			tmp.pos_b=o;
			tmp.ring=ring_index;
			tmp.index_in_ring=single_ring_num;
			n+=2;
			bump.push_back(tmp);
			cross_virtual_bump.push_back(tmp);
			single_ring_num++;
			}
		}
	}
	return ring_index;
}

void global_routing_sec1(vector<Driver> &bump, vector <vector <Driver> > &graph) {

	point tmp;
	char place[20];
	for (int i = 1; i<bump.size(); i++) {
		int current_ring = bump[i].ring;
		int current_ring_index = bump[i].index_in_ring - 1;
		int current_index = i - 1;
		int order = 1;
		int first = 0;
		int count = 0;
		int pos_a_last = bump[i].pos_a;
		int pos_b_last = bump[i].pos_b;

		tmp.index = bump[i].index;
		tmp.virtual_index = bump[i].virtual_index;
		while (current_ring>0 && bump[i].real) {
			while (current_ring_index>0) {
				if (order*(bump[i].virtual_index - bump[current_index].virtual_index)<0) {
					order = order*-1;
					if (i - current_index>1) {
						if (order == 1) {
							graph[bump[current_index].pos_a][bump[current_index].pos_b + 1].up.push_back(tmp);
							while ((pos_b_last - bump[current_index].pos_b)>0) {
								graph[bump[current_index].pos_a + 1][pos_b_last].left.push_back(tmp);
								pos_b_last -= 2;
							}
							pos_b_last = bump[current_index].pos_b;
						}
						else {
							graph[bump[current_index].pos_a][bump[current_index].pos_b + 1].down.push_back(tmp);
							while (pos_b_last - bump[current_index].pos_b>0) {
								graph[bump[current_index].pos_a - 1][pos_b_last].left.push_back(tmp);
								pos_b_last -= 2;
							}
							pos_b_last = bump[current_index].pos_b;
						}
					}
					first = 1;
				}
				current_index--;
				current_ring_index--;

			}

			if (order == -1) {
				graph[bump[current_index + 1].pos_a][bump[current_index + 1].pos_b - 1].up.push_back(tmp);
				while (pos_b_last - bump[current_index + 1].pos_b >= 0) {
					graph[bump[current_index + 1].pos_a + 1][pos_b_last].left.push_back(tmp);
					pos_b_last -= 2;
				}
				pos_b_last = bump[current_index + 1].pos_b;
				order = 1;
			}

			current_ring--;
			if (current_ring > 0) {
				current_ring_index = bump[current_index].index_in_ring;
				while (bump[current_index].virtual_index >  bump[i].virtual_index && current_ring_index>0) {
					current_index--;
					current_ring_index--;
				}

				if (current_ring_index != 0) {
					if (abs(bump[current_index].virtual_index - bump[i].virtual_index) >= 1)
						graph[bump[current_index].pos_a][bump[current_index].pos_b + 1].up.push_back(tmp);
					else {
						while (pos_b_last - bump[current_index].pos_b <= 0 && first == 0) {
							count = 1;
							graph[pos_a_last - 1][pos_b_last].right.push_back(tmp);
							pos_b_last += 2;
						}

						while (pos_b_last - bump[current_index].pos_b<0 && first == 1) {
							count = 1;
							graph[pos_a_last - 1][pos_b_last+2].right.push_back(tmp);
							pos_b_last += 2;
						}

						while (pos_b_last - bump[current_index].pos_b>=0 && count == 0 ) {
							graph[pos_a_last - 1][pos_b_last].left.push_back(tmp);
							pos_b_last -= 2;
						}

						break;
					}

					while (pos_b_last - bump[current_index].pos_b <= 0 && first == 0) {
						count = 1;
						graph[pos_a_last - 1][pos_b_last].right.push_back(tmp);
						pos_b_last += 2;
					}

					while (pos_b_last - bump[current_index].pos_b<0 && first == 1) {
						count = 1;
						graph[pos_a_last - 1][pos_b_last+2].right.push_back(tmp);
						pos_b_last += 2;
					}

					while (pos_b_last - bump[current_index].pos_b>0 && count == 0 ) {
						graph[pos_a_last - 1][pos_b_last].left.push_back(tmp);
						pos_b_last -= 2;
					}

					pos_a_last = bump[current_index].pos_a;
					pos_b_last = bump[current_index].pos_b;
				}

				else {
					graph[bump[current_index + 1].pos_a][bump[current_index + 1].pos_b - 1].up.push_back(tmp);
					while (pos_b_last - bump[current_index + 1].pos_b >= 0) {
						graph[pos_a_last - 1][pos_b_last].left.push_back(tmp);
						pos_b_last -= 2;
					}
					pos_a_last = bump[current_index + 1].pos_a;
					pos_b_last = bump[current_index + 1].pos_b - 2;
				}

				count = 0;
				first = 1;
			}
		}
	}
}

void global_routing_sec2(vector<Driver> &bump, vector <vector <Driver> > &graph) {

	point tmp;
	char place[20];
	for (int i = 1; i<bump.size(); i++) {
		int current_ring = bump[i].ring;
		int current_ring_index = bump[i].index_in_ring - 1;
		int current_index = i - 1;
		int order = 1;
		int count = 0;
		int first = 0;
		int pos_a_last = bump[i].pos_a;
		int pos_b_last = bump[i].pos_b;

		tmp.index = bump[i].index;
		tmp.virtual_index = bump[i].virtual_index;
		while (current_ring>0) {
			while (current_ring_index>0) {
				if (order*(bump[i].virtual_index - bump[current_index].virtual_index)<0) {
					order = order*-1;
					if (i - current_index>1) {
						if (order == 1) {
							graph[bump[current_index].pos_a - 1][bump[current_index].pos_b].right.push_back(tmp);
							while ((pos_a_last - bump[current_index].pos_a)>0) {
								graph[pos_a_last][bump[current_index].pos_b - 1].up.push_back(tmp);
								pos_a_last -= 2;
							}
							pos_a_last = bump[current_index].pos_a;
						}
						else {
							graph[bump[current_index].pos_a + 1][bump[current_index].pos_b].left.push_back(tmp);
							while (pos_a_last - bump[current_index].pos_a>0) {
								graph[pos_a_last][bump[current_index].pos_b + 1].up.push_back(tmp);
								pos_a_last -= 2;
							}
							pos_a_last = bump[current_index].pos_a;
						}
					}
					first = 1;
				}
				current_index--;
				current_ring_index--;

			}

			if (order == -1) {
				graph[bump[current_index + 1].pos_a - 1][bump[current_index + 1].pos_b].right.push_back(tmp);
				while (pos_a_last - bump[current_index + 1].pos_a >= 0) {
					graph[pos_a_last][bump[current_index].pos_b - 1].up.push_back(tmp);
					pos_a_last -= 2;
				}
				pos_a_last = bump[current_index + 1].pos_a;
				order = 1;
			}

			current_ring--;
			if (current_ring > 0) {
				current_ring_index = bump[current_index].index_in_ring;
				while (bump[current_index].virtual_index >  bump[i].virtual_index && current_ring_index>0) {
					current_index--;
					current_ring_index--;
				}

				if (current_ring_index != 0) {
					if (abs(bump[current_index].virtual_index - bump[i].virtual_index) >= 1)
						graph[bump[current_index].pos_a + 1][bump[current_index].pos_b].right.push_back(tmp);
					else {
						while (pos_a_last - bump[current_index].pos_a <= 0 && first == 0) {
							count = 1;
							graph[pos_a_last][pos_b_last + 1].down.push_back(tmp);
							pos_a_last += 2;
						}

						while (pos_a_last - bump[current_index].pos_a<0 && first == 1) {
							count = 1;
							graph[pos_a_last+2][pos_b_last + 1].down.push_back(tmp);
							pos_a_last += 2;
						}

						while (pos_a_last - bump[current_index].pos_a>=0 && count == 0) {
							graph[pos_a_last][pos_b_last + 1].up.push_back(tmp);
							pos_a_last -= 2;
						}
						break;
					}

					while (pos_a_last - bump[current_index].pos_a <= 0 && first == 0) {
						count = 1;
						graph[pos_a_last][pos_b_last + 1].down.push_back(tmp);
						pos_a_last += 2;
					}

					while (pos_a_last - bump[current_index].pos_a<0 && first == 1) {
						count = 1;
						graph[pos_a_last+2][pos_b_last + 1].down.push_back(tmp);
						pos_a_last += 2;
					}

					while (pos_a_last - bump[current_index].pos_a>0 && count == 0) {
						graph[pos_a_last][pos_b_last + 1].up.push_back(tmp);
						pos_a_last -= 2;
					}

					pos_a_last = bump[current_index].pos_a;
					pos_b_last = bump[current_index].pos_b;
				}

				else {
					graph[bump[current_index + 1].pos_a - 1][bump[current_index + 1].pos_b].right.push_back(tmp);
					while (pos_a_last - bump[current_index + 1].pos_a >= 0) {
						graph[pos_a_last][pos_b_last + 1].up.push_back(tmp);
						pos_a_last -= 2;
					}
					pos_a_last = bump[current_index + 1].pos_a - 2;
					pos_b_last = bump[current_index + 1].pos_b;

				}

				count = 0;
				first = 1;
			}
		}
	}
}

void global_routing_sec3(vector<Driver> &bump, vector <vector <Driver> > &graph) {

	point tmp;
	char place[20];
	for (int i = 1; i<bump.size(); i++) {
		int current_ring = bump[i].ring;
		int current_ring_index = bump[i].index_in_ring - 1;
		int current_index = i - 1;
		int order = 1;
		int count = 0;
		int first = 0;
		int pos_a_last = bump[i].pos_a;
		int pos_b_last = bump[i].pos_b;

		tmp.index = bump[i].index;
		tmp.virtual_index = bump[i].virtual_index;
		while (current_ring>0) {
			while (current_ring_index>0) {
				if (order*(bump[i].virtual_index - bump[current_index].virtual_index)<0) {
					order = order*-1;
					if (i - current_index>1) {
						if (order == 1) {
							graph[bump[current_index].pos_a][bump[current_index].pos_b - 1].down.push_back(tmp);
							while ((pos_b_last - bump[current_index].pos_b)<0) {
								graph[bump[current_index].pos_a - 1][pos_b_last].right.push_back(tmp);
								pos_b_last += 2;
							}
							pos_b_last = bump[current_index].pos_b;
						}
						else {
							graph[bump[current_index].pos_a][bump[current_index].pos_b - 1].up.push_back(tmp);
							while (pos_b_last - bump[current_index].pos_b<0) {
								graph[pos_a_last][bump[current_index].pos_b + 1].up.push_back(tmp);
								pos_b_last += 2;
							}
							pos_b_last = bump[current_index].pos_b;
						}
					}
					first = 1;
				}
				current_index--;
				current_ring_index--;

			}

			if (order == -1) {
				graph[bump[current_index + 1].pos_a][bump[current_index + 1].pos_b + 1].down.push_back(tmp);
				while (pos_b_last - bump[current_index + 1].pos_b <= 0) {
					graph[bump[current_index].pos_a - 1][pos_b_last].right.push_back(tmp);
					pos_b_last += 2;
				}
				pos_b_last = bump[current_index + 1].pos_b;
				order = 1;
			}

			current_ring--;
			if (current_ring > 0) {
				current_ring_index = bump[current_index].index_in_ring;
				while (bump[current_index].virtual_index >  bump[i].virtual_index && current_ring_index>0) {
					current_index--;
					current_ring_index--;
				}
				if (current_ring_index != 0) {
					if (abs(bump[current_index].virtual_index - bump[i].virtual_index) >= 1)
						graph[bump[current_index].pos_a][bump[current_index].pos_b - 1].down.push_back(tmp);
					else {
						while (pos_b_last - bump[current_index].pos_b >= 0 && first == 0) {
							count = 1;
							graph[pos_a_last + 1][pos_b_last].left.push_back(tmp);
							pos_b_last -= 2;
						}

						while (pos_b_last - bump[current_index].pos_b>0 && first == 1) {
							count = 1;
							graph[pos_a_last + 1][pos_b_last-2].left.push_back(tmp);
							pos_b_last -= 2;
						}

						while (pos_b_last - bump[current_index].pos_b<=0 && count == 0 ) {
							graph[pos_a_last + 1][pos_b_last].right.push_back(tmp);
							pos_b_last += 2;
						}
						break;
					}

					while (pos_b_last - bump[current_index].pos_b >= 0 && first == 0) {
						count = 1;
						graph[pos_a_last + 1][pos_b_last].left.push_back(tmp);
						pos_b_last -= 2;
					}

					while (pos_b_last - bump[current_index].pos_b>0 && first == 1) {
						count = 1;
						graph[pos_a_last + 1][pos_b_last-2].left.push_back(tmp);
						pos_b_last -= 2;
					}

					while (pos_b_last - bump[current_index].pos_b<0 && count == 0 ) {
						graph[pos_a_last + 1][pos_b_last].right.push_back(tmp);
						pos_b_last += 2;
					}

					pos_a_last = bump[current_index].pos_a;
					pos_b_last = bump[current_index].pos_b;
				}

				else {
					graph[bump[current_index + 1].pos_a][bump[current_index + 1].pos_b + 1].down.push_back(tmp);
					while (pos_b_last - bump[current_index + 1].pos_b <= 0) {
						graph[pos_a_last + 1][pos_b_last].right.push_back(tmp);
						pos_b_last += 2;
					}
					pos_a_last = bump[current_index + 1].pos_a;
					pos_b_last = bump[current_index + 1].pos_b + 2;
				}

				count = 0;
				first = 1;
			}
		}
	}
}

void global_routing_sec4(vector<Driver> &bump, vector <vector <Driver> > &graph) {

	point tmp;
	char place[20];
	for (int i = 1; i<bump.size(); i++) {
		int current_ring = bump[i].ring;
		int current_ring_index = bump[i].index_in_ring - 1;
		int current_index = i - 1;
		int first = 0;
		int order = 1;
		int count = 0;
		int pos_a_last = bump[i].pos_a;
		int pos_b_last = bump[i].pos_b;

		tmp.index = bump[i].index;
		tmp.virtual_index = bump[i].virtual_index;
		while (current_ring>0) {
			while (current_ring_index>0) {
				if (order*(bump[i].virtual_index - bump[current_index].virtual_index)<0) {
					order = order*-1;
					if (i - current_index>1) {
						if (order == 1) {
							graph[bump[current_index].pos_a + 1][bump[current_index].pos_b].left.push_back(tmp);
							while ((pos_a_last - bump[current_index].pos_a)<0) {
								graph[pos_a_last][bump[current_index].pos_b + 1].down.push_back(tmp);
								pos_a_last += 2;
							}
							pos_a_last = bump[current_index].pos_a;
						}
						else {
							graph[bump[current_index].pos_a - 1][bump[current_index].pos_b].right.push_back(tmp);
							while (pos_a_last - bump[current_index].pos_a<0) {
								graph[pos_a_last][bump[current_index].pos_b - 1].down.push_back(tmp);
								pos_a_last += 2;
							}
							pos_a_last = bump[current_index].pos_a;
						}
					}
					first = 1;
				}
				current_index--;
				current_ring_index--;

			}

			if (order == -1) {
				graph[bump[current_index + 1].pos_a + 1][bump[current_index + 1].pos_b].left.push_back(tmp);
				while (pos_a_last - bump[current_index + 1].pos_a <= 0) {
					graph[pos_a_last][bump[current_index].pos_b + 1].down.push_back(tmp);
					pos_a_last += 2;
				}
				pos_a_last = bump[current_index + 1].pos_a;
				order = 1;
			}

			current_ring--;
			if (current_ring > 0) {
				current_ring_index = bump[current_index].index_in_ring;
				while (bump[current_index].virtual_index >  bump[i].virtual_index && current_ring_index>0) {
					current_index--;
					current_ring_index--;
				}

				if (current_ring_index != 0) {
					if (abs(bump[current_index].virtual_index - bump[i].virtual_index) >= 1)
						graph[bump[current_index].pos_a - 1][bump[current_index].pos_b].left.push_back(tmp);
					else {
						while (pos_a_last - bump[current_index].pos_a >= 0 && first == 0) {
							count = 1;
							graph[pos_a_last][pos_b_last - 1].up.push_back(tmp);
							pos_a_last -= 2;
						}

						while (pos_a_last - bump[current_index].pos_a>0 && first == 1) {
							count = 1;
							graph[pos_a_last-2][pos_b_last - 1].up.push_back(tmp);
							pos_a_last -= 2;
						}


						while (pos_a_last - bump[current_index].pos_a<=0 && count == 0  ) {
							graph[pos_a_last][pos_b_last - 1].down.push_back(tmp);
							pos_a_last += 2;
						}
						break;
					}

					while (pos_a_last - bump[current_index].pos_a >= 0 && first == 0) {
						count = 1;
						graph[pos_a_last][pos_b_last - 1].up.push_back(tmp);
						pos_a_last -= 2;
					}

					while (pos_a_last - bump[current_index].pos_a>0 && first == 1) {
						count = 1;
						graph[pos_a_last-2][pos_b_last - 1].up.push_back(tmp);
						pos_a_last -= 2;
					}

					while (pos_a_last - bump[current_index].pos_a<0 && count == 0  ) {
						graph[pos_a_last][pos_b_last - 1].down.push_back(tmp);
						pos_a_last += 2;
					}

					pos_a_last = bump[current_index].pos_a;
					pos_b_last = bump[current_index].pos_b;
				}

				else {
					graph[bump[current_index + 1].pos_a + 1][bump[current_index + 1].pos_b].left.push_back(tmp);
					while (pos_a_last - bump[current_index + 1].pos_a <= 0) {
						graph[pos_a_last][pos_b_last - 1].down.push_back(tmp);
						pos_a_last += 2;
					}
					pos_a_last = bump[current_index + 1].pos_a + 2;
					pos_b_last = bump[current_index + 1].pos_b;
				}

				count = 0;
				first = 1;
			}
		}
	}
}

void global_routing_cross(vector<Driver> &bump, vector <vector <Driver> > &graph){

	point tmp;
	for(int i=0;i<bump.size();i++){
		tmp.index = bump[i].index;
		tmp.virtual_index = bump[i].virtual_index;
		for(int j=i+1;j<bump.size();j++){
			if(bump[i].index==bump[j].index){
				if(bump[i].pos_b-bump[j].pos_b>0){
					while(bump[i].pos_b-bump[j].pos_b!=0){
						graph[bump[i].pos_a+1][bump[i].pos_b].left.push_back(tmp);
						bump[i].pos_b-=2;
					}						
					while(bump[i].pos_a-bump[j].pos_a!=0){
						graph[bump[i].pos_a+2][bump[i].pos_b+1].down.push_back(tmp);
						bump[i].pos_a+=2;
					}
				}
				else if(bump[i].pos_b-bump[j].pos_b<0){
					while(bump[i].pos_b-bump[j].pos_b!=0){
						graph[bump[i].pos_a+1][bump[i].pos_b].right.push_back(tmp);
						bump[i].pos_b+=2;
					}
					while(bump[i].pos_a-bump[j].pos_a!=0){
						graph[bump[i].pos_a+2][bump[i].pos_b-1].down.push_back(tmp);
						bump[i].pos_a+=2;
					}
				}	
				else if(bump[i].pos_b-bump[j].pos_b==0){
					graph[bump[i].pos_a+1][bump[i].pos_b].right.push_back(tmp);
					while(bump[i].pos_a-bump[j].pos_a!=0){
						graph[bump[i].pos_a+2][bump[i].pos_b+1].down.push_back(tmp);
						bump[i].pos_a+=2;
					}
				}
			}
		}
	}
}

#endif
