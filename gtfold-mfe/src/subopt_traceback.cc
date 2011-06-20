/**
 GTfold: compute minimum free energy of RNA secondary structure
 Copyright (C) 2008  David A. Bader
 http://www.cc.gatech.edu/~bader

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 @author prashant {pgaurav@gatech.edu}

 */

#include <cstdio>

#include "constants.h"
#include "energy.h"
#include "utils.h"
#include "global.h"
#include "subopt_traceback.h"

const char* lstr[] = {"W", "V", "VBI", "VM", "WM"};

void (*trace_func[5]) (int i, int j, ps_t& ps, ps_stack_t& gs, int energy);

void print_stack(ps_stack_t temp)
{	
	std::cout << "{\n";
	while (!temp.empty())
	{
		ps_t ps = temp.top();
		temp.pop();
		ps.print();
		std::cout << '\n' ;
	}	
	std::cout << "}\n";
}

void process(ss_map_t& subopt_data, int energy, int len)
{
	int count = 0 ;
	ps_stack_t gstack;

	// initialize the partial structure
	// segment stack = {[1,n]}, label = W, list_bp = {} 
	ps_t first(0, len);
	first.push(segment(5, len, lW, W[len]));	
	gstack.push(first); // initialize the partial structure stacka

	while (1)
	{
		if (gstack.empty()) break; // exit

		ps_t ps = gstack.top();
		gstack.pop(); 

		if (ps.empty()) 
		{
			//std::cout << ps.str << " dG =" << ps.ae_ << std::endl;
			if (subopt_data.find(ps.str) == subopt_data.end())
			{
				subopt_data.insert(std::make_pair<std::string,int>(ps.str,ps.ae_));
				count++;
			}
			continue;
		}	
		else
		{
			segment smt = ps.top();
			ps.pop();
			size_t s1 = gstack.size();

			if (smt.j_ - smt.i_ >=4)
			{
				(*trace_func[smt.label_])(smt.i_, smt.j_, ps, gstack, energy);
			}

			if (gstack.size() == s1 && ps.total() == energy)
			{
				ps_t ps1(ps);
				gstack.push(ps1);
			}
		}
	}
}

ss_map_t subopt_traceback(int len, int delta)
{
	trace_func[0] = traceW;
	trace_func[1] = traceV;
	trace_func[2] = traceVBI;
	trace_func[3] = traceVM;
	//trace_func[4] = traceWM;

	int mfe = W[len];
	int level = 0; 

	ss_map_t subopt_data;

	for (level = mfe; level <= mfe + delta; level += 10)
		process(subopt_data, level, len);			

	return subopt_data;
}

void traceV(int i, int j, ps_t& ps, ps_stack_t& gstack, int energy)
{
	// Hairpin Loop
	if (eH(i,j) + ps.total()  == energy )
	{
		//	std::cout << "Hairpin " << i  << ' ' << j << std::endl;
		ps_t ps1(ps); 
		ps1.accumulate(eH(i,j));
		ps1.update(i, j, '(', ')');
		push_to_gstack(gstack, ps1);
	}

	// Stack
	if (eS(i, j) + V(i+1, j-1) + ps.total() == energy)
	{
		//	std::cout << "Stack " << i  << ' ' << j << std::endl;
		ps_t ps1(ps);
		ps1.push(segment(i+1, j-1, lV, V(i+1, j-1)));
		ps1.accumulate(eS(i,j));
		ps1.update(i, j , '(', ')');
		push_to_gstack(gstack, ps1);
	}

	// Internal Loop
	if (VBI(i,j) + ps.total() == energy )
	{
		//std::cout << "Internal " << i  << ' ' << j << std::endl;
		ps_t ps1(ps);
		ps1.push(segment(i, j, lVBI, VBI(i,j)));
		ps1.update(i, j, '(', ')');
		push_to_gstack(gstack, ps1);
	}

	// Multiloop
	if ( VM(i,j) + ps.total() == energy )
	{
		//	std::cout << "Multi " << i  << ' ' << j << std::endl;
		ps_t ps1(ps);
		ps1.push(segment(i, j, lVM, VM(i,j)));
		ps1.update(i, j, '(', ')');
		push_to_gstack(gstack, ps1);
	}

}

void traceVBI(int i, int j, ps_t& ps, ps_stack_t& gstack, int energy)
{
	int p,q;
	for (p = i+1; p < MIN(j-2-TURN, i+MAXLOOP+1); ++p)
	{
		int minq = j-i+p-MAXLOOP-2;
		if (minq < p+1+TURN) minq = p+1+TURN;

		for (q = minq; q < j; q++) 
		{
			if (V(p, q) + eL(i, j, p, q) + ps.total() == energy )
			{
				ps_t ps1(ps);
				ps1.push(segment(p, q, lV, V(p, q)));
				ps1.accumulate(eL(i, j, p, q));
				push_to_gstack(gstack, ps1);
			}
		}
	}
}

void traceW(int h, int j, ps_t& ps, ps_stack_t& gstack, int energy)
{
	if ( (h>=j) || j == 0 || j == 1) return;

	for (int i = h; i < j; ++i)
	{
		int wim1 =  MIN(0, W[i-1]);

		int wij = V(i,j) + auPenalty(i, j) + wim1;
		if (wij + ps.total() == energy )
		{
			ps_t ps1(ps);
			ps1.push(segment(i, j, lV, V(i,j)));
			if (wim1 <= 0 && i > 1) ps1.push(segment(1, i-1, lW, W[i-1]));
			ps1.accumulate(auPenalty(i, j));
			push_to_gstack(gstack, ps1);
		}		

		int wijd = V(i,j-1) + auPenalty(i,j-1) + Ed5(j-1,i,j) + wim1;
		if (wijd + ps.total() == energy )
		{
			ps_t ps3(ps);
			ps3.push(segment(i, j-1, lV, V(i, j-1)));
			if (wim1 <= 0 && i > 1)	ps3.push(segment(1, i-1, lW, W[i-1]));
			ps3.accumulate(auPenalty(i,j-1) + Ed5(j-1,i,j)); 
			push_to_gstack(gstack, ps3);
		}	

		int widj = V(i+1,j) + auPenalty(i+1,j) + Ed3(j,i+1,i) + wim1;
		if (widj + ps.total() == energy )
		{
			ps_t ps4(ps);
			ps4.push(segment(i+1, j, lV, V(i+1,j)));
			if (wim1 <= 0 && i > 1) ps4.push(segment(1, i-1, lW, W[i-1]));
			ps4.accumulate(auPenalty(i+1, j) + Ed3(j,i+1,i));
			push_to_gstack(gstack, ps4);
		}

		int widjd = V(i+1, j-1) + auPenalty(i+1,j-1) + Ed3(j-1, i+1, i) + Ed5(j-1, i+1, j) + wim1;
		if (widjd + ps.total() == energy )
		{
			ps_t ps2(ps);
			ps2.push(segment(i+1, j-1, lV, V(i+1,j-1)));
			if (wim1 <= 0 && i > 1) ps2.push(segment(1, i-1, lW, W[i-1]));
			ps2.accumulate(auPenalty(i+1,j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j));
			push_to_gstack(gstack, ps2);
		}
	}

	if (W[j-1] + ps.total() == energy )
	{
		ps_t ps1(ps);
		ps1.push(segment(1, j-1, lW, W[j-1]));
		push_to_gstack(gstack, ps1);
	}	
}

void traceWM(ps_t& ps, ps_map_t& filter, int energy)
{
	ps_stack_t wm_stack;
	wm_stack.push(ps);

	while (!wm_stack.empty())
	{
		ps_t pss = wm_stack.top();
		wm_stack.pop();

		// check if all segments in pss are Vs
		if (pss.empty()) 
		{
			std::string key = pss.str;
			while (!pss.empty_v())
			{
				segment ss = pss.top_v();
				pss.pop_v();

				pss.push(ss);
				key[ss.i_-1] = '*'; key[ss.j_-1] = '*';
			}
			filter.insert(std::pair<std::string, ps_t>(key, pss));
			continue;
		}

		segment wm_seg = pss.top();
		pss.pop();

		int i1 = wm_seg.i_; int j1 = wm_seg.j_;

		if (wm_seg.label_ == lV)
		{
			ps_t ps1(pss);
			ps1.push_v(wm_seg);
			wm_stack.push(ps1);
			continue;
		}
		int h = i1; int k = j1;

		int wmij = V(h,k) + auPenalty(h,k) + Eb;
		if (pss.total() + wmij == energy )
		{
			ps_t ps1(pss);
			ps1.push_v(segment(h,k,lV, V(h,k)));
			ps1.accumulate(auPenalty(h,k) + Eb);
			wm_stack.push(ps1);
		}

		int wmijd = V(h,k-1) + Ed5(k-1,h,k)+ auPenalty(h,k-1) +Eb+ Ec ;
		if (pss.total() + wmijd == energy )
		{
			ps_t ps1(pss);
			ps1.push_v(segment(h,k-1,lV, V(h,k-1)));
			ps1.accumulate(Ed5(k-1,h,k)+ auPenalty(h,k-1) +Eb+Ec);
			wm_stack.push(ps1);
		}

		int wmidj = V(h+1,k) + Ed3(k,h+1,h) + auPenalty(h+1, k) + Eb+Ec ;
		if (pss.total() + wmidj == energy )
		{
			ps_t ps1(pss);
			ps1.push_v(segment(h+1, k, lV, V(h+1,k)));
			ps1.accumulate(Ed3(k,h+1,h) + auPenalty(h+1,k) +Eb+Ec);
			wm_stack.push(ps1);
		}

		int wmidjd = V(h+1, k-1) + Ed3(k-1,h+1,h) + Ed5(k-1,h+1,k) +  auPenalty(h+1,k-1) +Eb+ 2*Ec ;
		if (pss.total() + wmidjd == energy )
		{
			ps_t ps1(pss);
			ps1.push_v(segment(h+1, k-1, lV, V(h+1,k-1)));
			ps1.accumulate(Ed3(k-1,h+1,h) + Ed5(k-1,h+1,k) + auPenalty(h+1,k-1) +Eb+ 2*Ec);
			wm_stack.push(ps1);
		}


		if (pss.total() + WM(i1,j1-1)  + Ec  == energy )
		{
			ps_t ps1(pss);
			ps1.push(segment(i1,j1-1, lWM, WM(i1,j1-1)));		
			ps1.accumulate(Ec);
			wm_stack.push(ps1);
		}

		if (pss.total() + WM[i1+1][j1] + Ec == energy )
		{
			ps_t ps1(pss);
			ps1.push(segment(i1+1,j1, lWM, WM(i1+1,j1)));		
			ps1.accumulate(Ec);
			wm_stack.push(ps1);
		}

		for (int h = i1+1; h <= j1-1; ++h)
		{	
			if (WM(i1,h) + WM(h+1,j1) + pss.total() == energy )
			{
				ps_t ps1(pss);
				ps1.push(segment(i1, h, lWM, WM(i1,h)));
				ps1.push(segment(h+1, j1, lWM, WM(h+1,j1)));
				wm_stack.push(ps1);
			}
		}
	}
}

void traceVM(int i, int j, ps_t& ps, ps_stack_t& gstack, int energy)
{
	std::map<std::string, ps_t> filter;
	int h;	

	for (h = i+1; h <= j-1; ++h)
	{
		int dG;
		int d5 = Ed5(i,j,i+1);
		int d3 = Ed3(i,j,j-1);
		int common = auPenalty(i,j) + Ea + Eb;

		dG = common + WM[i+1][h-1] + WM[h][j-1];
		if (dG + ps.total()  == energy )
		{
			ps_t ps1(ps);
			ps1.push(segment(i+1,h-1, lWM, WM[i+1][h-1]));
			ps1.push(segment(h, j-1, lWM, WM[h][j-1]));
			ps1.accumulate(common) ; 
			traceWM(ps1, filter, energy);
		}

		dG =  common + WM[i+2][h-1] + WM[h][j-1] + d5 + Ec;
		if (dG + ps.total()  == energy )
		{
			ps_t ps1(ps);
			ps1.push(segment(i+2,h-1, lWM, WM[i+2][h-1]));	
			ps1.push(segment(h, j-1, lWM, WM[h][j-1]));	
			ps1.accumulate(common + d5 + Ec);
			traceWM(ps1, filter, energy);
		}

		dG = common + WM[i+1][h-1] + WM[h][j-2] + d3 + Ec;
		if (dG + ps.total()  == energy )
		{
			ps_t ps1(ps);	
			ps1.push(segment(i+1,h-1, lWM, WM[i+1][h-1]));	
			ps1.push(segment(h, j-2, lWM, WM[h][j-2]));	
			ps1.accumulate(common + d3 + Ec);
			traceWM(ps1, filter, energy);
		}

		dG = common + WM[i+2][h-1] + WM[h][j-2] + d5 + d3 + 2*Ec;
		if (dG + ps.total()  == energy )
		{ 	
			ps_t ps1(ps);	
			ps1.push(segment(i+2,h-1, lWM, WM[i+2][h-1]));	
			ps1.push(segment(h, j-2, lWM, WM[h][j-2]));	
			ps1.accumulate(common + d3 + d5 + 2*Ec) ;
			traceWM(ps1, filter, energy);
		}
	}

	ps_map_t::iterator it ;
	for (it = filter.begin(); it != filter.end(); ++it)
	{
		push_to_gstack(gstack, it->second);
	}
}

void push_to_gstack(ps_stack_t& gstack , const ps_t& v)
{
	gstack.push(v);
}
