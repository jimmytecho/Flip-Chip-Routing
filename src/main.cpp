#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "class.h"
#include "function.h"

using namespace std;

int main(int argc, char **argv){

	ifstream fin(argv[1]);
	if(argc!=2) {
		cout<<"Wrong format!!"<<endl;
		return 0;
	}
	if (!fin.is_open()){
		cout<<"Cannot open input file "<<argv[1]<<"!!"<<endl;
		return 0;
	}
	
 	//read global data
		int bx0,by0,bxf,byf; //boundary of the entire design
		int width,spacing,driver_num,bump_num;
		string junk;
		fin >> junk >> bx0 >> by0 >> bxf >> byf;
		fin >> junk >> width;
		fin >> junk >> spacing;
		fin >> junk >> driver_num;
		vector <Driver> driver(driver_num);
 	//END read global data

	//read detailed data
		for (int i=0; i< driver_num ; i++)
			fin >> driver[i].index >> driver[i].x0 >> driver[i].y0 >> driver[i].xf >> driver[i].yf;
		fin >> junk >> bump_num;
		vector <Driver> bump(bump_num);
		for (int i=0; i< bump_num ; i++)
			fin >> bump[i].index >> bump[i].x0 >> bump[i].y0 >>bump[i].xf >> bump[i].yf;
		fin.close();
	//END read detailed data

	//assign virtual_index of each bump, initialize with index
		for(int i=0;i<bump.size();i++){
			bump[i].virtual_index=bump[i].index;
		}
	//END assign virtual_index of each bump, initialize with index

	//open vector <vector <Driver> > &graph
		int x0_min=1000000000,y0_min=1000000,x0_max=0,y0_max=0;
		for(int i=0;i<bump_num;i++){
			if(bump[i].x0<x0_min || bump[i].y0<y0_min){
				x0_min=bump[i].x0;
				y0_min=bump[i].y0;
			}
			if(bump[i].x0>x0_max || bump[i].y0>y0_max){
				y0_max=bump[i].y0;
				x0_max=bump[i].x0;
			}
		}

		int x1_min=10000000,y1_min=100000000,x1_max=0,y1_max=0;
		for(int i=0;i<bump_num;i++){
			if((bump[i].x0<x1_min || bump[i].y0<y1_min)&& bump[i].x0!=x0_min){
				x1_min=bump[i].x0;
				y1_min=bump[i].y0;
			}
			if((bump[i].x0>x1_max || bump[i].y0>y1_max)&& bump[i].y0!=y0_max){
				y1_max=bump[i].y0;
				x1_max=bump[i].x0;
			}
		}

		int side1=(x0_max-x0_min)/(x1_min-x0_min);
		int side2=(y0_max-y0_min)/(y0_max-y1_max);
		int xdis=(x1_min-x0_min);
		int ydis=(y0_max-y1_max);
		int wi=2*(side1+1)+1;
		int le=2*(side2+1)+1;

		//open a two-dimensional vector: graph 
		vector<vector<Driver> > graph (wi, vector<Driver> (le));

		int xlength=bump[1].xf-bump[1].x0;//x邊長
		int ylength=bump[1].yf-bump[1].y0;//y邊長

		for(int i=0;i<le;i=i+2){
			for(int j=0;j<wi;j=j+2){
				graph[i][j].x0=(x0_min-xdis+xlength)+(j/2)*xdis;
				graph[i][j].y0=(y0_max+ylength)-(i/2)*(ydis);
				graph[i][j].xf=(x0_min)+(j/2)*(xdis);
				graph[i][j].yf=(y0_max+ydis)-(i/2)*(ydis);
			}
		}

		for(int j=0;j<le;j++){
			for(int i=(j+1)%2;i<wi;i=i+2){
				graph[j][i].pos_a=j;
				graph[j][i].pos_b=i;
			}
		}

		x0_min-=xdis;
		y0_max+=ydis;

		for(int i=0;i<bump_num;i++){  
			int w= (bump[i].x0-x0_min)/xdis;
			int l= (y0_max-bump[i].y0)/ydis;

			bump[i].pos_a=2*l-1;
			bump[i].pos_b=2*w-1;
		}
	//END open vector <vector <Driver> > &graph

	//bump[i].index  <->  i
		vector <int> indexd(driver_num+1);
		for (int i = 1; i <= driver_num; i++) indexd[driver[i-1].index] = i-1;
		vector <int> indexb(bump_num);
		for (int j = 1; j <= bump_num; j++) indexb[bump[j-1].index] = j-1;
	//END bump[i].index  <->  i

	//assign driver's section
		int countl=0,countr=0,countu=0,countd=0;
		for(int i=0; i<driver_num; i++){
			if (driver[i].yf == byf){
				driver[i].sec.push_back(1);
				driver[i].sec.push_back(1);
				countu++;
			}
			else if (driver[i].x0 == bx0){
				driver[i].sec.push_back(4);
				driver[i].sec.push_back(4);
				countl++;
			}
			else if (driver[i].xf == bxf){
				driver[i].sec.push_back(2);
				driver[i].sec.push_back(2);
				countr++;
			}
			else if (driver[i].y0 == by0){
				driver[i].sec.push_back(3);
				driver[i].sec.push_back(3);
				countd++;
			}
			bump[i].sec.push_back(0);
			bump[i].sec.push_back(0);
		}
	//END assign driver's section
	
	//partition  set bump[]sec[0],sec[1]
		vector <Driver> crossed_line_bump1;//a vector to store the bumps on the left_up crossed line
		vector <Driver> crossed_line_bump2;//a vector to store the bumps on the right_up crossed line
		vector <Driver> crossed_line_bump3;//a vector to store the bumps on the left_down crossed line
		vector <Driver> crossed_line_bump4;//a vector to store the bumps on the right_down crossed line
		double midx,midy;
		midx = ((double) bx0 + (double) bxf)/2.0;//the middlex of the graph
		midy = ((double) by0 + (double) byf)/2.0;//the middley of the graph
		for (int i=0; i<bump_num; i++){
			double newyf = (double) bump[i].yf;
			double newx0 = (double) bump[i].x0; 
			double newxf = (double) bump[i].xf;
			double newy0 = (double) bump[i].y0;
			double middlex, middley, m1, m2, mp, mn;
			middlex = (newx0+newxf)/2.0;//middle of the bump
			middley = (newy0+newyf)/2.0;//middle of the bump
			mp = (byf-middley)/(bxf-middlex);//positive slope of the graph
			mn = (byf-middley)/(bx0-middlex);//negative slope of the graph
			if (middley >= midy){
				if (middlex < midx){
					m1 = (newy0-midy)/(newx0-midx);//the left slope of the bump
					m2 = (newyf-midy)/(newx0-midx);//the right slope of the bump
					if (m1<mn && m2<mn ) bump[i].sec[0] = 1;//both slopes are on the left of mn
					else if (m1>mn && m2>mn) bump[i].sec[0] = 4;//both slopes are on the right of mn
					else {//means the bump is on the crossed line
						crossed_line_bump1.push_back(bump[i]);
						if (driver[indexd[bump[i].index]].sec[0] == 1) bump[i].sec[0] = 1;//if the driver is in section 1, then set bump to section 1
						else if (driver[indexd[bump[i].index]].sec[0] == 4) bump[i].sec[0] = 4;//if the driver is in section 4, then set bump to section 4
					}
				}     
				else if (middlex > midx){
					m1 = (newy0-midy)/(newxf-midx);
					m2 = (newyf-midy)/(newx0-midx);
					if (m1>mp && m2>mp) bump[i].sec[0] = 1;
					else if (m1<mp && m2<mp) bump[i].sec[0] = 2;
					else {
						crossed_line_bump2.push_back(bump[i]);
						if (driver[indexd[bump[i].index]].sec[0] == 1) bump[i].sec[0] = 1;
						else if (driver[indexd[bump[i].index]].sec[0] == 2) bump[i].sec[0] = 2;
					}
				}
				else if (middlex == midx) bump[i].sec[0] = 1;//if the bump is on the vertical line
			}

			else if (middley < midy){
				if (middlex < midx){
					m1 = (newyf-midy)/(newx0-midx);
					m2 = (newy0-midy)/(newxf-midx);
					if (m1>mp && m2>mp) bump[i].sec[0] = 3;
					else if (m1<mp && m2<mp) bump[i].sec[0] = 4;
					else {
						crossed_line_bump3.push_back(bump[i]);            
						if (driver[indexd[bump[i].index]].sec[0] == 3) bump[i].sec[0] = 3;
						else if (driver[indexd[bump[i].index]].sec[0] == 4) bump[i].sec[0] = 4;
					}
				}
				else if (middlex > midx){
					m1 = (newyf-midy)/(newx0-midx);
					m2 = (newy0-midy)/(newxf-midx);
					if (m1<mn && m2<mn) bump[i].sec[0] = 3;
					else if (m1>mn && m2>mn) bump[i].sec[0] = 2;
					else {
						crossed_line_bump4.push_back(bump[i]);
						if (driver[indexd[bump[i].index]].sec[0] == 3) bump[i].sec[0] = 3;
						else if (driver[indexd[bump[i].index]].sec[0] == 2) bump[i].sec[0] = 2;
					}
				}
				else if (middlex == midx) bump[i].sec[0] = 3;
			} 
		} 
	//END partition  set bump[]sec[0],sec[1]

	//find out the cross section bump
		int cross_num = 0;
		vector <Driver> pos_line;//store original ross-section bumps
		for (int i=0;i<bump_num;i++){
			if (bump[i].sec[0] != driver[indexd[bump[i].index]].sec[0]){
				bump[i].sec[1] = driver[indexd[bump[i].index]].sec[0];
				driver[indexd[bump[i].index]-1].sec[1] = bump[i].sec[0];
				pos_line.push_back(bump[i]);
				cross_num++;
			}
			else bump[i].sec[1] = bump[i].sec[0];
		}
	//END find out the cross section bump
  
	//set graph[i][j].index & graph[i][j].sec[0] & graph[i][j].sec[1]
		for (int i=1;i<le;i=i+2){
			for (int j=1;j<wi;j=j+2){
				graph[i][j].sec.push_back(0);
				graph[i][j].sec.push_back(0);
			}
		}
		for (int i=0;i<bump.size();i++){
			graph[bump[i].pos_a][bump[i].pos_b].index = bump[i].index;
			graph[bump[i].pos_a][bump[i].pos_b].virtual_index = bump[i].virtual_index;
			graph[bump[i].pos_a][bump[i].pos_b].sec[0] = bump[i].sec[0];
			graph[bump[i].pos_a][bump[i].pos_b].sec[1] = bump[i].sec[1];
		}
	//END set graph[i][j].index & graph[i][j].sec[0] & graph[i][j].sec[1]

	//if cross section can be reduced to sectional
		for(int i=0;i<pos_line.size();i++){
			if(pos_line[i].sec[0]==1 && pos_line[i].sec[1]==2){
				int a=graph[crossed_line_bump2[i].pos_a][crossed_line_bump2[i].pos_b+2].sec[1];
				int b=graph[crossed_line_bump2[i].pos_a+2][crossed_line_bump2[i].pos_b].sec[1];
				if (a==2 && b==2) bump[indexb[pos_line[i].index]].sec[0]=2;
			}
			else if(pos_line[i].sec[0]==1 && pos_line[i].sec[1]==4){
				int a=graph[crossed_line_bump2[i].pos_a][crossed_line_bump2[i].pos_b-2].sec[1];
				int b=graph[crossed_line_bump2[i].pos_a+2][crossed_line_bump2[i].pos_b].sec[1];
				if (a==4 && b==4) bump[indexb[pos_line[i].index]].sec[0]=4;
			}
			else if(pos_line[i].sec[0]==2 && pos_line[i].sec[1]==1){
				int a=graph[crossed_line_bump2[i].pos_a][crossed_line_bump2[i].pos_b-2].sec[1];
				int b=graph[crossed_line_bump2[i].pos_a-2][crossed_line_bump2[i].pos_b].sec[1];
				if (a==1 && b==1) bump[indexb[pos_line[i].index]].sec[0]=1;
			}
			else if(pos_line[i].sec[0]==2 && pos_line[i].sec[1]==3){
				int a=graph[crossed_line_bump2[i].pos_a][crossed_line_bump2[i].pos_b-2].sec[1];
				int b=graph[crossed_line_bump2[i].pos_a+2][crossed_line_bump2[i].pos_b].sec[1];
				if (a==3 && b==3) bump[indexb[pos_line[i].index]].sec[0]=3;
			}
			else if(pos_line[i].sec[0]==3 && pos_line[i].sec[1]==4){
				int a=graph[crossed_line_bump2[i].pos_a][crossed_line_bump2[i].pos_b-2].sec[1];
				int b=graph[crossed_line_bump2[i].pos_a-2][crossed_line_bump2[i].pos_b].sec[1];
				if (a==4 && b==4) bump[indexb[pos_line[i].index]].sec[0]=4;
			}
			else if(pos_line[i].sec[0]==3 && pos_line[i].sec[1]==2){
				int a=graph[crossed_line_bump2[i].pos_a][crossed_line_bump2[i].pos_b+2].sec[1];
				int b=graph[crossed_line_bump2[i].pos_a-2][crossed_line_bump2[i].pos_b].sec[1];
				if (a==2 && b==2) bump[indexb[pos_line[i].index]].sec[0]=2;
			}
			else if(pos_line[i].sec[0]==4 && pos_line[i].sec[1]==3){
				int a=graph[crossed_line_bump2[i].pos_a][crossed_line_bump2[i].pos_b+2].sec[1];
				int b=graph[crossed_line_bump2[i].pos_a+2][crossed_line_bump2[i].pos_b].sec[1];
				if (a==3 && b==3) bump[indexb[pos_line[i].index]].sec[0]=3;
			}
			else if(pos_line[i].sec[0]==4 && pos_line[i].sec[1]==1){
				int a=graph[crossed_line_bump2[i].pos_a][crossed_line_bump2[i].pos_b+2].sec[1];
				int b=graph[crossed_line_bump2[i].pos_a-2][crossed_line_bump2[i].pos_b].sec[1];
				if (a==1 && b==1) bump[indexb[pos_line[i].index]].sec[0]=1;
			}
		}
	//END if cross section can be reduced to sectional

	//create four vectors to do global routing
		vector<Driver> bumper_sec1;
		vector<Driver> bumper_sec2;
		vector<Driver> bumper_sec3;
		vector<Driver> bumper_sec4;

		for(int i=0;i<bump.size();i++){
			if(bump[i].sec[0]==1) bumper_sec1.push_back(bump[i]);
			if(bump[i].sec[0]==2) bumper_sec2.push_back(bump[i]);
			if(bump[i].sec[0]==3) bumper_sec3.push_back(bump[i]);
			if(bump[i].sec[0]==4) bumper_sec4.push_back(bump[i]);
		}

		int size1=bumper_sec1.size();
		int size2=bumper_sec2.size();
		int size3=bumper_sec3.size();
		int size4=bumper_sec4.size();
		vector<Driver> cross_virtual_bump;
	//END create four vectors to do global routing

	//setring
		int max_ring1 = 0;
		int max_ring2 = 0;
		int max_ring3 = 0;
		int max_ring4 = 0;

		Sort_x(bumper_sec1,"O");
		Sort_y(bumper_sec1,"U");
		max_ring1 = SetRing(bumper_sec1,"U",cross_virtual_bump);

		Sort_y(bumper_sec2,"U");
		Sort_x(bumper_sec2,"U");
		max_ring2 = SetRing(bumper_sec2,"R",cross_virtual_bump);

		Sort_x(bumper_sec3,"U");
		Sort_y(bumper_sec3,"O");
		max_ring3 = SetRing(bumper_sec3,"D",cross_virtual_bump);

		Sort_y(bumper_sec4,"O");
		Sort_x(bumper_sec4,"O");
		max_ring4 = SetRing(bumper_sec4,"L",cross_virtual_bump);
	//END setring

	//virtual bump in destination section
		for(int i=0;i<size1;i++){
			if(bumper_sec1[i].sec[0]!=bumper_sec1[i].sec[1]){
				if(bumper_sec1[i].sec[1]==2){
					Driver cross_bumper;
					if(bumper_sec2.back().ring!=max_ring1){
						cross_bumper.index=bumper_sec1[i].index;
						cross_bumper.ring=max_ring1;
						cross_bumper.index_in_ring=1;
						cross_bumper.pos_a = bumper_sec2.back().pos_a-10;
						cross_bumper.pos_b = bumper_sec2.back().pos_b-2;
						cross_bumper.virtual_index=bumper_sec1[i].index;
						cross_bumper.sec.push_back(bumper_sec1[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec1[i].sec[1]);
					}
					else{
						cross_bumper.index=bumper_sec1[i].index;
						cross_bumper.ring=max_ring1;
						cross_bumper.index_in_ring=bumper_sec2.back().index_in_ring+1;
						cross_bumper.pos_a = bumper_sec2.back().pos_a+2;
						cross_bumper.pos_b = bumper_sec2.back().pos_b;
						cross_bumper.virtual_index=bumper_sec1[i].index;
						cross_bumper.sec.push_back(bumper_sec1[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec1[i].sec[1]);
					}
					bumper_sec2.push_back(cross_bumper);
					cross_virtual_bump.push_back(cross_bumper);
				}

				if(bumper_sec1[i].sec[1]==3){
					Driver cross_bumper;
					if(bumper_sec3.back().ring!=max_ring1){
						cross_bumper.index=bumper_sec1[i].index;
						cross_bumper.ring=max_ring1;
						cross_bumper.index_in_ring=1;
						cross_bumper.pos_a = bumper_sec3.back().pos_a-2;
						cross_bumper.pos_b = bumper_sec3.back().pos_b+8;
						cross_bumper.virtual_index=bumper_sec1[i].index;
						cross_bumper.sec.push_back(bumper_sec1[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec1[i].sec[1]);
					}
					else{
						cross_bumper.index=bumper_sec1[i].index;
						cross_bumper.ring=max_ring1;
						cross_bumper.index_in_ring=bumper_sec3.back().index_in_ring+1;
						cross_bumper.pos_a = bumper_sec3.back().pos_a;
						cross_bumper.pos_b = bumper_sec3.back().pos_b+2;
						cross_bumper.virtual_index=bumper_sec1[i].index;
						cross_bumper.sec.push_back(bumper_sec1[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec1[i].sec[1]);
					}
					bumper_sec3.push_back(cross_bumper);
					cross_virtual_bump.push_back(cross_bumper);   
				}

				if(bumper_sec1[i].sec[1]==4){
					Driver cross_bumper;
					if(bumper_sec4.back().ring!=max_ring1){
						cross_bumper.index=bumper_sec1[i].index;
						cross_bumper.ring=max_ring1;
						cross_bumper.index_in_ring=1;
						cross_bumper.pos_a = bumper_sec4.back().pos_a+10;
						cross_bumper.pos_b = bumper_sec4.back().pos_b+2;
						cross_bumper.virtual_index=bumper_sec1[i].index;
						cross_bumper.sec.push_back(bumper_sec1[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec1[i].sec[1]);
					}
					else{
						cross_bumper.index=bumper_sec1[i].index;
						cross_bumper.ring=max_ring1;
						cross_bumper.index_in_ring=bumper_sec4.back().index_in_ring+1;
						cross_bumper.pos_a = bumper_sec4.back().pos_a+2;
						cross_bumper.pos_b = bumper_sec4.back().pos_b;
						cross_bumper.virtual_index=bumper_sec1[i].index;
						cross_bumper.sec.push_back(bumper_sec1[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec1[i].sec[1]);
						}
					bumper_sec4.push_back(cross_bumper);
					cross_virtual_bump.push_back(cross_bumper);
				}
			}
		}

		for(int i=0;i<size2;i++){
			if(bumper_sec2[i].sec[0]!=bumper_sec2[i].sec[1]){
				if(bumper_sec2[i].sec[1]==3){
					Driver cross_bumper;
					if(bumper_sec3.back().ring!=max_ring2){
						cross_bumper.index=bumper_sec2[i].index;
						cross_bumper.ring=max_ring2;
						cross_bumper.index_in_ring=1;
						cross_bumper.pos_a = bumper_sec3.back().pos_a-2;
						cross_bumper.pos_b = bumper_sec3.back().pos_b+8;
						cross_bumper.virtual_index=bumper_sec2[i].index;
						cross_bumper.sec.push_back(bumper_sec2[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec2[i].sec[1]);
					}
					else{
						cross_bumper.index=bumper_sec2[i].index;
						cross_bumper.ring=max_ring2;
						cross_bumper.index_in_ring=bumper_sec3.back().index_in_ring+1;
						cross_bumper.pos_a = bumper_sec3.back().pos_a;
						cross_bumper.pos_b = bumper_sec3.back().pos_b+2;
						cross_bumper.virtual_index=bumper_sec2[i].index;
						cross_bumper.sec.push_back(bumper_sec2[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec2[i].sec[1]);
					}
					bumper_sec3.push_back(cross_bumper);
					cross_virtual_bump.push_back(cross_bumper); 
				}

				if(bumper_sec2[i].sec[1]==4){
					Driver cross_bumper;
					if(bumper_sec4.back().ring!=max_ring2){
						cross_bumper.index=bumper_sec2[i].index;
						cross_bumper.ring=max_ring2;
						cross_bumper.index_in_ring=1;
						cross_bumper.pos_a = bumper_sec3.back().pos_a+10;
						cross_bumper.pos_b = bumper_sec3.back().pos_b+2;
						cross_bumper.virtual_index=bumper_sec2[i].index;
						cross_bumper.sec.push_back(bumper_sec2[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec2[i].sec[1]);
					}
					else{
						cross_bumper.index=bumper_sec2[i].index;
						cross_bumper.ring=max_ring2;
						cross_bumper.index_in_ring=bumper_sec4.back().index_in_ring+1;
						cross_bumper.pos_a = bumper_sec3.back().pos_a+2;
						cross_bumper.pos_b = bumper_sec3.back().pos_b;
						cross_bumper.virtual_index=bumper_sec2[i].index;
						cross_bumper.sec.push_back(bumper_sec2[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec2[i].sec[1]);
					}
					bumper_sec4.push_back(cross_bumper);
					cross_virtual_bump.push_back(cross_bumper);
				}

				if(bumper_sec2[i].sec[1]==1){
					Driver cross_bumper;
					if(bumper_sec1.back().ring!=max_ring2){
						cross_bumper.index=bumper_sec2[i].index;
						cross_bumper.ring=max_ring2;
						cross_bumper.index_in_ring=1;
						cross_bumper.pos_a = bumper_sec1.back().pos_a+2;
						cross_bumper.pos_b = bumper_sec1.back().pos_b-8;
						cross_bumper.virtual_index=bumper_sec2[i].index;
						cross_bumper.sec.push_back(bumper_sec2[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec2[i].sec[1]);
					}
					else{
						cross_bumper.index=bumper_sec2[i].index;
						cross_bumper.ring=max_ring2;
						cross_bumper.index_in_ring=bumper_sec1.back().index_in_ring+1;
						cross_bumper.pos_a = bumper_sec1.back().pos_a;
						cross_bumper.pos_b = bumper_sec1.back().pos_b+2;
						cross_bumper.virtual_index=bumper_sec2[i].index;
						cross_bumper.sec.push_back(bumper_sec2[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec2[i].sec[1]);
					}
					bumper_sec1.push_back(cross_bumper);
					cross_virtual_bump.push_back(cross_bumper);
				}
			}
		}

		for(int i=0;i<size3;i++){
			if(bumper_sec3[i].sec[0]!=bumper_sec3[i].sec[1]){
				if(bumper_sec3[i].sec[1]==4){
					Driver cross_bumper;
					if(bumper_sec4.back().ring!=max_ring3){
						cross_bumper.index=bumper_sec3[i].index;
						cross_bumper.ring=max_ring3;
						cross_bumper.index_in_ring=1;
						cross_bumper.pos_a = bumper_sec4.back().pos_a+10;
						cross_bumper.pos_b = bumper_sec4.back().pos_b+2;
						cross_bumper.virtual_index=bumper_sec3[i].index;
						cross_bumper.sec.push_back(bumper_sec3[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec3[i].sec[1]);
					}
					else{
						cross_bumper.index=bumper_sec3[i].index;
						cross_bumper.ring=max_ring3;
						cross_bumper.index_in_ring=bumper_sec4.back().index_in_ring+1;
						cross_bumper.pos_a = bumper_sec4.back().pos_a+2;
						cross_bumper.pos_b = bumper_sec4.back().pos_b;
						cross_bumper.virtual_index=bumper_sec3[i].index;
						cross_bumper.sec.push_back(bumper_sec3[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec3[i].sec[1]);
					}
					bumper_sec4.push_back(cross_bumper);
					cross_virtual_bump.push_back(cross_bumper);
				}

				if(bumper_sec3[i].sec[1]==1){
					Driver cross_bumper;
					if(bumper_sec1.back().ring!=max_ring3){
						cross_bumper.index=bumper_sec3[i].index;
						cross_bumper.ring=max_ring3;
						cross_bumper.index_in_ring=1;
						cross_bumper.pos_a = bumper_sec1.back().pos_a+2;
						cross_bumper.pos_b = bumper_sec1.back().pos_b-8;
						cross_bumper.virtual_index=bumper_sec3[i].index;
						cross_bumper.sec.push_back(bumper_sec3[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec3[i].sec[1]);
					}
					else{
						cross_bumper.index=bumper_sec3[i].index;
						cross_bumper.ring=max_ring3;
						cross_bumper.index_in_ring=bumper_sec1.back().index_in_ring+1;
						cross_bumper.pos_a = bumper_sec1.back().pos_a;
						cross_bumper.pos_b = bumper_sec1.back().pos_b+2;
						cross_bumper.virtual_index=bumper_sec3[i].index;
						cross_bumper.sec.push_back(bumper_sec3[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec3[i].sec[1]);
					}
					bumper_sec1.push_back(cross_bumper);
					cross_virtual_bump.push_back(cross_bumper);
				}

				if(bumper_sec3[i].sec[1]==2){
					Driver cross_bumper;
					if(bumper_sec2.back().ring!=max_ring3){
						cross_bumper.index=bumper_sec3[i].index;
						cross_bumper.ring=max_ring3;
						cross_bumper.index_in_ring=1;
						cross_bumper.pos_a = bumper_sec2.back().pos_a-10;
						cross_bumper.pos_b = bumper_sec2.back().pos_b-2;
						cross_bumper.virtual_index=bumper_sec3[i].index;
						cross_bumper.sec.push_back(bumper_sec3[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec3[i].sec[1]);
					}
					else{
						cross_bumper.index=bumper_sec3[i].index;
						cross_bumper.ring=max_ring3;
						cross_bumper.index_in_ring=bumper_sec2.back().index_in_ring+1;
						cross_bumper.pos_a = bumper_sec2.back().pos_a+2;
						cross_bumper.pos_b = bumper_sec2.back().pos_b;
						cross_bumper.virtual_index=bumper_sec3[i].index;
						cross_bumper.sec.push_back(bumper_sec3[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec3[i].sec[1]);
					}
					bumper_sec2.push_back(cross_bumper);
					cross_virtual_bump.push_back(cross_bumper);
				}
			}
		}

		for(int i=0;i<size4;i++){
			if(bumper_sec4[i].sec[0]!=bumper_sec4[i].sec[1]){
				if(bumper_sec4[i].sec[1]==1){
					Driver cross_bumper;
					if(bumper_sec1.back().ring!=max_ring4){
						cross_bumper.index=bumper_sec4[i].index;
						cross_bumper.ring=max_ring4;
						cross_bumper.index_in_ring=1;
						cross_bumper.pos_a = bumper_sec1.back().pos_a+2;
						cross_bumper.pos_b = bumper_sec1.back().pos_b-8;
						cross_bumper.virtual_index=bumper_sec4[i].index;
						cross_bumper.sec.push_back(bumper_sec4[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec4[i].sec[1]);
					}
					else{
						cross_bumper.index=bumper_sec4[i].index;
						cross_bumper.ring=max_ring4;
						cross_bumper.index_in_ring=bumper_sec1.back().index_in_ring+1;
						cross_bumper.pos_a = bumper_sec1.back().pos_a;
						cross_bumper.pos_b = bumper_sec1.back().pos_b+2;
						cross_bumper.virtual_index=bumper_sec4[i].index;
						cross_bumper.sec.push_back(bumper_sec4[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec4[i].sec[1]);
					}
					bumper_sec1.push_back(cross_bumper);
					cross_virtual_bump.push_back(cross_bumper);
				}

				if(bumper_sec4[i].sec[1]==2){
					Driver cross_bumper;
					if(bumper_sec2.back().ring!=max_ring4){
						cross_bumper.index=bumper_sec4[i].index;
						cross_bumper.ring=max_ring4;
						cross_bumper.index_in_ring=1;
						cross_bumper.pos_a = bumper_sec2.back().pos_a-2;
						cross_bumper.pos_b = bumper_sec2.back().pos_b-10;
						cross_bumper.virtual_index=bumper_sec4[i].index;
						cross_bumper.sec.push_back(bumper_sec4[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec4[i].sec[1]);
					}
					else{
						cross_bumper.index=bumper_sec4[i].index;
						cross_bumper.ring=max_ring4;
						cross_bumper.index_in_ring=bumper_sec2.back().index_in_ring+1;
						cross_bumper.pos_a = bumper_sec2.back().pos_a+2;
						cross_bumper.pos_b = bumper_sec2.back().pos_b;
						cross_bumper.virtual_index=bumper_sec4[i].index;
						cross_bumper.sec.push_back(bumper_sec4[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec4[i].sec[1]);
					}
					bumper_sec2.push_back(cross_bumper);
					cross_virtual_bump.push_back(cross_bumper);
				}

				if(bumper_sec4[i].sec[1]==3){
					Driver cross_bumper;
					if(bumper_sec3.back().ring!=max_ring4){
						cross_bumper.index=bumper_sec4[i].index;
						cross_bumper.ring=max_ring4;
						cross_bumper.index_in_ring=1;
						cross_bumper.pos_a = bumper_sec3.back().pos_a-2;
						cross_bumper.pos_b = bumper_sec3.back().pos_b+8;
						cross_bumper.virtual_index=bumper_sec4[i].index;
						cross_bumper.sec.push_back(bumper_sec4[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec4[i].sec[1]);
					}
					else{
						cross_bumper.index=bumper_sec4[i].index;
						cross_bumper.ring=max_ring4;
						cross_bumper.index_in_ring=bumper_sec3.back().index_in_ring+1;
						cross_bumper.pos_a = bumper_sec3.back().pos_a;
						cross_bumper.pos_b = bumper_sec3.back().pos_b-2;
						cross_bumper.virtual_index=bumper_sec4[i].index;
						cross_bumper.sec.push_back(bumper_sec4[i].sec[1]);
						cross_bumper.sec.push_back(bumper_sec4[i].sec[1]);
					}
					bumper_sec3.push_back(cross_bumper);
					cross_virtual_bump.push_back(cross_bumper);
				}
			}
		}
	//END virtual bump in destination section

	//store the original cross_virtual_bump before it's changed by other function
		vector <Driver> new_virtual_bump;
		new_virtual_bump = cross_virtual_bump;

	//global routing, cross-section too
		global_routing_sec1(bumper_sec1,graph);
		global_routing_sec2(bumper_sec2,graph);
		global_routing_sec3(bumper_sec3,graph);
		global_routing_sec4(bumper_sec4,graph);
		global_routing_cross(cross_virtual_bump,graph);
	//END endglobal routing, cross-section too

	//set graph[i][j] for virtual bump
		for (int i=0;i<bumper_sec1.size();i++){
			graph[bumper_sec1[i].pos_a][bumper_sec1[i].pos_b].index = bumper_sec1[i].index;
			graph[bumper_sec1[i].pos_a][bumper_sec1[i].pos_b].virtual_index = bumper_sec1[i].virtual_index;
		}
		for (int i=0;i<bumper_sec2.size();i++){
			graph[bumper_sec2[i].pos_a][bumper_sec2[i].pos_b].index = bumper_sec2[i].index;
			graph[bumper_sec2[i].pos_a][bumper_sec2[i].pos_b].virtual_index = bumper_sec2[i].virtual_index;
		}
		for (int i=0;i<bumper_sec3.size();i++){
			graph[bumper_sec3[i].pos_a][bumper_sec3[i].pos_b].index = bumper_sec3[i].index;
			graph[bumper_sec3[i].pos_a][bumper_sec3[i].pos_b].virtual_index = bumper_sec3[i].virtual_index;
		}
		for (int i=0;i<bumper_sec4.size();i++){
			graph[bumper_sec4[i].pos_a][bumper_sec4[i].pos_b].index = bumper_sec4[i].index;
			graph[bumper_sec4[i].pos_a][bumper_sec4[i].pos_b].virtual_index = bumper_sec4[i].virtual_index;
		}
	//END set graph[i][j] for virtual bump

	//fill
		for (int i=2;i<le-2;i=i+2) {     
			for (int j=2;j<wi-2;j=j+2) fill(graph,i,j);
		}
	//END fill

	//detailed routing for bump interval
		vector <path> mypath(bump_num+1);
		for(int i=2;i<le-2;i=i+2){
			for(int j=2;j<wi-2;j=j+2) detailed(graph,mypath,i,j,width,spacing);
		}
		for (int i=0;i<new_virtual_bump.size();i++){
			int xx0 = x0_min+xdis;
			int xxf = x0_min+xdis+xlength;
			int yy0 = y0_max-ydis;
			int yyf = y0_max-ydis+ylength;

			line newt;
			newt.x1 = xx0 + (new_virtual_bump[i].pos_b-1)*xdis/2;
			newt.y1 = yy0 - (new_virtual_bump[i].pos_a-1)*ydis/2;
			newt.x2 = xxf + (new_virtual_bump[i].pos_b-1)*xdis/2;
			newt.y2 = yy0 - (new_virtual_bump[i].pos_a-1)*ydis/2;
			mypath[new_virtual_bump[i].index].path.push_back(newt);
			line newt2;
			newt2.x1 = xxf + (new_virtual_bump[i].pos_b-1)*xdis/2;
			newt2.y1 = yy0 - (new_virtual_bump[i].pos_a-1)*ydis/2;
			newt2.x2 = xxf + (new_virtual_bump[i].pos_b-1)*xdis/2;
			newt2.y2 = yyf - (new_virtual_bump[i].pos_a-1)*ydis/2;
			mypath[new_virtual_bump[i].index].path.push_back(newt2);
		}
	//END detailed routing for bump interval

	//linking for tile
		for (int i = 3; i < le-2; i = i + 2) {
			for (int j = 2; j < wi-2; j = j + 2) link(graph, mypath, i, j, width, spacing);
		}

		for (int i = 2; i < le-2; i = i + 2) {
			for (int j = 3; j < wi-2; j = j + 2) link2(graph, mypath, i, j, width, spacing);
		}
	//END linking for tile

	//connecting driver to the most outer part
		vector<Driver> driver_sec1;
		vector<Driver> driver_sec2;
		vector<Driver> driver_sec3;
		vector<Driver> driver_sec4;

		for (int i = 0; i<driver.size(); i++) {
			if (driver[i].sec[0] == 1) driver_sec1.push_back(driver[i]);
			if (driver[i].sec[0] == 2) driver_sec2.push_back(driver[i]);
			if (driver[i].sec[0] == 3) driver_sec3.push_back(driver[i]);
			if (driver[i].sec[0] == 4) driver_sec4.push_back(driver[i]);
		}

		for (int i = 0; i < driver_sec1.size(); i++) {
			indexd[driver_sec1[i].index] = i;
			indexd[driver_sec2[i].index] = i;
			indexd[driver_sec3[i].index] = i;
			indexd[driver_sec4[i].index] = i;
		}
		int AAA;
		for (int i = 0;i<bump.size();i++){
			AAA = bump[i].index;
			indexb[AAA] = i;
		}

		initialize(driver_sec1, bump, indexd, indexb, 1, wi, le);
		Walk(graph, driver_sec1, bump, indexd, indexb, mypath, width, spacing);
		initialize(driver_sec2, bump, indexd, indexb, 2, wi, le);
		Walk2(graph, driver_sec2, bump, indexd, indexb, mypath, width, spacing, wi);
		initialize(driver_sec3, bump, indexd, indexb, 3, wi, le);
		Walk3(graph, driver_sec3, bump, indexd, indexb, mypath, width, spacing, le);
		initialize(driver_sec4, bump, indexd, indexb, 4, wi, le);
		Walk4(graph, driver_sec4, bump, indexd, indexb, mypath, width, spacing);
	//END connecting driver to the most outer part

	//linking for the outer part;
		for (int i = 0; i < mypath[144].path.size(); i++) {
			line L;
			L = mypath[144].path[i];
		}

		for (int j = 2; j < wi - 2; j = j + 2) link(graph, mypath, 1, j, width, spacing);
		for (int i = 2; i < le - 2; i = i + 2) link2(graph, mypath, i, wi - 2, width, spacing);
		for (int j = 2; j < wi - 2; j = j + 2) link(graph, mypath, le - 2, j, width, spacing);
		for (int i = 2; i < le - 2; i = i + 2) link2(graph, mypath, i, 1, width, spacing);
	//END linking for the outer part;
   
	//fout
		ofstream fout;
		fout.open("../outputs/outputfile.txt");

		for(int i=1;i<=144;i++){
			int aaa=i;
			for (int i = 0; i<mypath[aaa].path.size(); i++) {
				fout <<aaa<<" "<< mypath[aaa].path[i].x1 << " " << mypath[aaa].path[i].y1<<" "<<mypath[aaa].path[i].x2 << " " << mypath[aaa].path[i].y2 << endl;
			}
		}

		for (int i = 0; i < driver.size(); i++) {
			fout << i << " " << driver[i].x0 << " " << driver[i].y0 << " " << driver[i].x0 << " " << driver[i].yf << endl;
			fout << i << " " << driver[i].x0 << " " << driver[i].y0 << " " << driver[i].xf << " " << driver[i].y0 << endl;
			fout << i << " " << driver[i].xf << " " << driver[i].y0 << " " << driver[i].xf << " " << driver[i].yf << endl;
			fout << i << " " << driver[i].x0 << " " << driver[i].yf << " " << driver[i].xf << " " << driver[i].yf << endl;
		}
		for (int i = 0; i < bump.size(); i++) {
			fout << i << " " << bump[i].x0 << " " << bump[i].y0 << " " << bump[i].x0 << " " << bump[i].yf << endl;
			fout << i << " " << bump[i].x0 << " " << bump[i].y0 << " " << bump[i].xf << " " << bump[i].y0 << endl;
			fout << i << " " << bump[i].xf << " " << bump[i].y0 << " " << bump[i].xf << " " << bump[i].yf << endl;
			fout << i << " " << bump[i].x0 << " " << bump[i].yf << " " << bump[i].xf << " " << bump[i].yf << endl;
		}

		fout.close();
	//END fout
	return 0;
    
}
