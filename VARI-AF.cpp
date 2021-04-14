#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <libgen.h> // basename
#include <cassert>
#include "tclap/CmdLine.h"
#include <map>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#include "debruijn_graph_shifted.hpp"
#include "algorithm.hpp"
#include "cosmo-color.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
//#include <regex.h>
//using namespace std;
//REF_LEN/using namespace sdsl;

#include <sys/timeb.h>
#include <cmath>
#include <unistd.h>
//#include "colored-dbg.hpp"
//#pragma comment(linker, "/STACK:1024000000")
//struct rlimit lim;
//if(getrlimit(RLIMIT_STACK, &lim) < 0){ 
//perror("getrlimit"); 
//exit(1); 
//} 
//printf("stack size:%lu\n",lim.rlim_cur);
//static char base[] = {'?','A','C','G','T'};

static char base[] = {'$','A','C','G','T'};
//static kmer = dbg.k - 1;

//static map<ssize_t,int> nodepos;
static debruijn_graph_shifted<> dbg;
static sd_vector<> colors;
//static vector<int> nestlevelset;
static vector<int> mycolors;
//static int REF_LEN = 16569;
static string researched_samples;
static string dbgstartstring = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGG";
static string dbgendstring = "TCTGGTTCCTACTTCAGGGTCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATG";
//static string dbgstartstring = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA";
vector<int> plusnodes;
//static vector<int> mycolors(colors.size(),0);
//static vector<int> mycoverage(colors.size(),0);
//static vector<int> myplusstrand(dbg.num_edges(),0);

int getMilliCount(){
  timeb tb;
  ftime(&tb);
  int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  return nCount;
}


int getMilliSpan(int nTimeStart){
  int nSpan = getMilliCount() - nTimeStart;
  if(nSpan < 0)
    nSpan += 0x100000 * 1000;
  return nSpan;
}

void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd(BANNER, ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input", ".dbg file.", true, "", "graph_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> color_filename_arg("color", ".rrr file.", true, "", "color_file", cmd);
  string color_mask1 = "color_mask1";
  TCLAP::ValueArg<std::string> color_mask1_arg("a", "color_mask1",
	    "Color mask 1, color1 [" + color_mask1 + "]", false, "", color_mask1, cmd);
  string color_mask2 = "color_mask2";
  TCLAP::ValueArg<std::string> color_mask2_arg("b", "color_mask2",
	    "Color mask 2, color2 [" + color_mask2 + "]", false, "", color_mask2, cmd);
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.color_filename  = color_filename_arg.getValue();
  params.color_mask1     = color_mask1_arg.getValue();
  params.color_mask2     = color_mask2_arg.getValue();
}

void test_symmetry(debruijn_graph_shifted<> dbg) {
  for (unsigned long x = 0; x<dbg.sigma+1;x++) {
    ssize_t in = dbg.incoming(43, x);
    if (in == -1)
      continue;
    for (unsigned long y = 0; y<dbg.sigma+1;y++) {
      ssize_t out = dbg.outgoing(in, y);
      if (out == -1)
	continue;
      cout << "Incoming " << in <<  ":" << out <<"\n";
    }
  }
}

template <typename T>
T Max(vector<T>& vec){
	if(vec.size() == 0)
		return -1;
	for(size_t i = 1;i < vec.size();i++){
		if(vec[0] < vec[i]){
	    	vec[0] = vec[i];
		}
	}
	return vec[0];
}

template <typename T>
T Min(vector<T>& vec){
	if(vec.size() == 0)
		return -1;
	for(size_t i = 1;i < vec.size();i++){
		if(vec[0] > vec[i]){
	    	vec[0] = vec[i];
		}
	}
	return vec[0];
}


struct Basis{
	int type;//1:supernode;2:branches
	vector<ssize_t> incomingnode;
	vector<ssize_t> outgoingnode;
	ssize_t startnode;
	ssize_t endnode;
	string color;
	int length;
	int width;
	string sequence;
};

struct BubInfo{
	int bubbleid,level,maxlevel,numnodes,numedges,numbranches,numsupernodes;
	string bubcolor;
	vector<ssize_t> pairnodeset;//(level,nodes1,nodes2,...)
	vector<int> affiliation;//(bub_id,father_id,child_id1,child_id2)
	vector<int> colorcount;
	vector<string> branches;
	vector<vector<ssize_t>> branchnodes;
	vector<vector<string>> branchcolors;
	vector<string> bridgebranches,tipbranches;
	vector<ssize_t> sinknodes,sourcenodes;
	vector<vector<int>> branchdigits;
	vector<vector<ssize_t>> bridgebranchnodes,tipbranchnodes;
//	ssize_t firstnode;//Cyclic iteration's startnode
	vector<int> visited;
	vector<Basis> bubblebasis;
	vector<vector<int>> weightmat;
	vector<int> sampleweight,samplelength;
	int maxlength;
	vector<int> gapbasisid;
	string variantspos;
	vector<string> pathseqs,pathcolors;
	vector<int> pathlengths,pathweights,pathsupernodes;
	vector<vector<int>> pathbranchwidths;
	map<ssize_t,int> nodespos,refnodepos;
};

int StringCMP(string str1,string str2){//-1:diff length; 0:same length && no common 1 value; 1:same length && have common 1 value but less than length; 2:same
	if(str1 == str2)
		return 2;
	int str1len = str1.size();
	int str2len = str2.size();
	if(str1len != str2len)
		return -1;
	int samecount = 0;
	for(int i = 0;i < str1len;i++){
		if(str1[i] == str2[i] && str1[i]!='0')
			samecount++;
	}
	if(samecount == 0)
		return 0;
	if(samecount > 0 && samecount < str1len)
		return 1;
	else{
		cout << "string comparing meets error:\n" << str1 << "\n" << str2 << endl;
		return -1;
	}
}
//#####################################################################################################
// To get the color of the node
string Getcolor(ssize_t node){
	int numcolors = colors.size()/dbg.size();
	string zerostr(numcolors,'0');
	string colstr = zerostr;
//	cout << node << endl;
	for(unsigned long x = 1;x < 5;x++){
		string tmpcolstr;
		ssize_t nextnode = dbg.outgoing(node,x);
		if(nextnode == -1)
			continue;
		ssize_t nextedge = dbg.outgoing_edge(node,x);
		bool findcoming = false;
		for(unsigned long x2 = 0;x2 < 5;x2++){
			ssize_t frontnode = dbg.incoming(nextnode,x2);
			if(frontnode == -1)
				continue;
			if(frontnode == node){
				findcoming = true;
				nextedge = dbg.incoming_edge(nextnode,x2);
				break;
			}
		}
		if(findcoming == false)
			cout << "doesn't find the front node" << endl;
		for(int c = 0; c < numcolors; c++){
			stringstream ss;
			int tmpint = colors[nextedge * numcolors + c];
			ss << tmpint;
			string tmpstr = ss.str();
			tmpcolstr += tmpstr;
		}
//		if(zerostr == tmpcolstr){
//			string subcolor = Getcolor(nextnode);
//			tmpcolstr = subcolor;
//		}
		int len1 = colstr.size();
		int len2 = tmpcolstr.size();
		if(len1 == numcolors && len2 == numcolors && colstr != tmpcolstr){
			for(int i = 0;i < numcolors;i++){
				if(!(colstr[i] == '0' && tmpcolstr[i] == '0')){
					colstr[i] = '1';
				}
			}
		}
		else{
			colstr = tmpcolstr;
		}
	}
//	cout << node << "'s color:" << colstr << endl;
	return colstr;
}

string getnodecolor(ssize_t node){
	int numcolors = colors.size()/dbg.size();
//	cout << node <<endl;
	string zerostr(numcolors,'0');
	string colstr = zerostr;
	for(unsigned long x = 0;x < 5;x++){//get the nodes' color from the incoming edges
		string tmpcolstr;
		ssize_t edge;
		edge = dbg.incoming_edge(node,x);
		if (edge == -1)
			continue;
//		cout << edge << endl;
		for(int c = 0; c < numcolors; c++){
			stringstream ss;
			int tmpint = colors[edge * numcolors + c];
			ss << tmpint;
			string tmpstr = ss.str();
			tmpcolstr += tmpstr;
		}
		int len1 = colstr.size();
		int len2 = tmpcolstr.size();
//		cout << len1 << "\t" << len2 << "\t" << numcolors << "\t" << colstr << "\t" << tmpcolstr << endl;
		if(len1 == numcolors && len2 == numcolors && colstr != tmpcolstr){
			for(int i = 0;i < numcolors;i++){
				if(!(colstr[i] == '0' && tmpcolstr[i] == '0')){
					colstr[i] = '1';
				}
			}
		}
		else{
			colstr = tmpcolstr;
		}
//		cout << "getnodecolor:" << node << "->" << colstr << endl;
//		if(zerostr == colstr){
//			colstr = Getcolor(node);
//		}
	}
	if(dbg.indegree(node) == 0)
		colstr = Getcolor(node);
	return colstr;
}
//#####################################################################################################
// To get the color of the edge
string getedgecolor(ssize_t edge){
	string colstr;
	int numcolors = colors.size()/dbg.size();
	for(int c = 0; c < numcolors; c++){
		stringstream ss;
		int tmpint = colors[edge * numcolors + c];
		ss << tmpint;
		string tmpstr = ss.str();
		colstr += tmpstr;
	}
	return colstr;
}
//#####################################################################################################
// To get the new degree of the node
int getindegree(vector<int> nestlevelset,int nestinglevel,ssize_t node){
	int indegree = 0;
	for(unsigned long i = 0;i < dbg.sigma + 1;i++){
		ssize_t edge = dbg.incoming_edge(node,i);
		if(edge == -1)
			continue;
		ssize_t frontnode = dbg.incoming(node,i);
		if(nestlevelset[frontnode] == nestinglevel - 1 && nestinglevel > 1)
			continue;
		indegree ++;
	}
	return indegree;
}

int getoutdegree(vector<int> nestlevelset,int nestinglevel,ssize_t node){
	int outdegree = 0;
	for(unsigned long i = 1;i < dbg.sigma + 1;i++){
		ssize_t edge = dbg.outgoing_edge(node,i);
		if(edge == -1)
			continue;
		ssize_t nextnode = dbg.outgoing(node,i);
		if(nestlevelset[nextnode] == nestinglevel - 1 && nestinglevel > 1)
			continue;
		outdegree ++;
	}
	return outdegree;
}
//#####################################################################################################
//#####################################################################################################
//#####################################################################################################
// To get the max or min level of the nestlevelset
int GetMaxLevel(vector<int> nestlevelset){
	int max = 1;
	ssize_t num_nodes = dbg.num_nodes();
	for(ssize_t node = 0;node < num_nodes;node++){
		if(nestlevelset[node] > max)
			max = nestlevelset[node];
	}
	return max;
}
int GetMinLevel(vector<int> nestlevelset){
	int min = 20;
	ssize_t num_nodes = dbg.num_nodes();
	for(ssize_t node = 0;node < num_nodes;node++){
		if(nestlevelset[node] < min && nestlevelset[node] > 0)
			min = nestlevelset[node];
	}
	return min;
}

//#####################################################################################################
string ColorCap(vector<string> colorset){
	int num = colorset.size();
	if(num < 1){
		cout << "color set length is error"  << endl;
		return "0";
	}
	string capstring = colorset[0];
	int len = colorset[0].size();
	for(int i = 1;i < num;i++){
		string str = colorset[i];
		int tmplen = str.size();
		if(tmplen == len){
			for(int j = 0;j < len;j++){
				if(str[j] == '1' && capstring[j] == '1')
					capstring[j] = '1';
				else{
					capstring[j] = '0';
				}
			}
		}
		else{
			cout << "string sets' length is not the same" << endl;
			return "0";
		}
	}
	return capstring;
}

string ColorCup(vector<string> colorset){
	int num = colorset.size();
//	cout << "There are " << num << " strings:";
	int numcolors = colors.size()/dbg.size();
	string cupstring(numcolors,'0');
	if(num < 1){
		cout << "color set length is error"  << endl;
		return cupstring;
	}
	cupstring = colorset[0];
	int len = colorset[0].size();
	for(int i = 1;i < num;i++){
		string str = colorset[i];
		int tmplen = str.size();
		if(tmplen == len){
			for(int j = 0;j < len;j++){
				if(str[j] == '1' || cupstring[j] == '1')
					cupstring[j] = '1';
			}
		}
		else{
			cout << "string sets' length is not the same" << endl;
			return cupstring;
		}
	}
	return cupstring;
}

string VarCup(vector<string> varset){
	int num = varset.size();
//	cout << "There are " << num << " strings:";
//	int numcolors = colors.size()/dbg.size();
//	string cupstring(numcolors,'0');
	string cupstring = "";
	if(num < 1){
		cout << "variants set length is error"  << endl;
		return cupstring;
	}
//	cupstring = varset[0];
//	int len = colorset[0].size();
	int len = 0;
	for(int i = 0;i < num;i++){
		string str = varset[i];
		if(str == "")
			continue;
		if(cupstring == "")
			cupstring = str;
		int tmplen = str.size();
		if(len == 0 || tmplen == len){
			for(int j = 0;j < len;j++){
				if(str[j] == '1' || cupstring[j] == '1')
					cupstring[j] = '1';
				if(str[j] == '2' || cupstring[j] == '2')
					cupstring[j] = '2';
				if((str[j] == '1' && cupstring[j] == '2') || (str[j] == '2' && cupstring[j] == '1')){
					cupstring[j] = '3';
					
				}
			}
			if(len == 0)
				len = tmplen;
		}
		else{
			cout << "variants string sets' length is not the same" << endl;
			return cupstring;
		}
	}
	return cupstring;
}

int myindegree(ssize_t node,string color){
	int indegree = 0;
	int numcolors = colors.size()/dbg.size();
	string zerocolor(numcolors,'0');
	for(int k = 1;k < 5;k++){
		ssize_t frontnode = dbg.incoming(node,k);
		if(frontnode == -1)
			continue;
		string frontnodecolor = getnodecolor(frontnode);
		string capcolor = ColorCap({researched_samples,color,frontnodecolor});
		if(capcolor == zerocolor)
			continue;
		indegree++;		
	}
	return indegree;
}
int myoutdegree(ssize_t node,string color){
	int outdegree = 0;
	int numcolors = colors.size()/dbg.size();
	string zerocolor(numcolors,'0');
	for(int k = 1;k < 5;k++){
		ssize_t nextnode = dbg.outgoing(node,k);
		if(nextnode == -1)
			continue;
		string nextnodecolor = getnodecolor(nextnode);
		string capcolor = ColorCap({researched_samples,color,nextnodecolor});
		if(capcolor == zerocolor)
			continue;
		outdegree++;		
	}
	return outdegree;
}

const char *const starts[] = {"GCCATACTGCGTCATGTCGCCCTGACGCGC","GCAGGTTCGAATCCTGCACGACCCACCAAT","GCTTAACCTCACAACCCGAAGATGTTTCTT","AAAACCCGCCGAAGCGGGTTTTTACGTAAA","AATCCTGCACGACCCACCAGTTTTAACATC","AGAGTTCCCCGCGCCAGCGGGGATAAACCG","GAATACGTGCGCAACAACCGTCTTCCGGAG"};

//#####################################################################################################################################
//##################################################  SuperBubbles Analysis  ##########################################################
//#####################################################################################################################################
vector<int> GetPlusnode(){
	int t = getMilliCount();
	int kmer = dbg.k-1;
	int count = 0;
	string startlabel = dbgstartstring.substr(1,kmer);
	//cout << startlabel << endl;
	ssize_t numnodes = dbg.num_nodes();
	vector<int> plusnode(numnodes,0);
	ssize_t startnode = 0;
	int numcolors = colors.size()/dbg.size();
	string stdcolor(numcolors,'1');
//	cout << stdcolor << endl;
	for(ssize_t node = 0;node < numnodes;node++){
		string tmpnodelabel = dbg.node_label(node);
		if(tmpnodelabel == startlabel){
//			cout << startlabel << endl;
			string nodecol = getnodecolor(node);
//			cout << nodecol << endl;
			if(nodecol == stdcolor){
				startnode = node;
				cout << "find start node: " << node << endl;
				break;
			}
		}
	}
	queue<ssize_t> nodelist;
	nodelist.push(startnode);
	plusnode[startnode] = 1;
	while(!nodelist.empty()){
		ssize_t node = nodelist.front();
		count++;
//		cout << node << endl;
		for(int k = 1;k < 5;k++){
			ssize_t nextnode = dbg.outgoing(node,k);
			if(nextnode == -1||plusnode[nextnode] == 1)
				continue;
			if(plusnode[nextnode] == 1){cout << "There are some problems while getting plus node" << endl;}
			while(dbg.indegree(nextnode) == 1 && dbg.outdegree(nextnode) == 1){
				plusnode[nextnode] = 1;
				count++;
				for(int k1 = 1;k1 < 5;k1++){
					ssize_t pos = dbg.outgoing(nextnode,k1);
					if(pos == -1)
						continue;
					nextnode = pos;
					break;
				}
			//	cout << " " << nextnode << endl;
			}
			if(plusnode[nextnode] == 0){
				plusnode[nextnode] = 1;
				nodelist.push(nextnode);
			}
		}
		nodelist.pop();
	}

	cout << "The number of plus node: " << count << endl; 
	int totaltime = getMilliSpan(t);
	if(totaltime < 1000){
		cout << "Total calculation time: " << getMilliSpan(t) << "ms" << endl;
	}
	else{
		int tolsecond = totaltime/1000;
		int hour = tolsecond/3600;
		int second1 = tolsecond%3600;
		int minute = second1/60;
		int second = second1%60;
		cout << "Total get strand information's time: " << hour << "h" << minute << "m" << second << "s" << endl;
	}
	cerr << "Finished getting strand information" << endl;
	return plusnode;
}

//vector<vector<int>> Get_Cut_Pos(map<ssize_t,int> nodepos,vector<map<ssize_t,int>>samp_nodepos,vector<vector<ssize_t>> cycle_visiting_nodes){///输出每个环在每个样本中的切割位点
void Get_Cut_Pos(map<ssize_t,int> nodepos,vector<vector<ssize_t>> cycle_visiting_nodes){
//vector<vector<int>> cycle_samp_cutpos = Get_Cut_Pos(map<ssize_t,int> nodepos,vector<map<ssize_t,int>>samp_nodepos,vector<vector<ssize_t>> cycle_visiting_nodes){///输出每个环在每个样本中的切割位点
	//int cycle_num = sorted_cycle_visiting_nodes.size();
	int cycle_num = cycle_visiting_nodes.size();
	ofstream cutfile("mycutinfo.txt");
	if(cycle_num > 0){//需要判断是否存在交集
		cout << "There are " << cycle_num << " cycles" << endl;
//		如果各自无交，那么不做任何处理，在各自的区域内寻找切割点；
//		如果有包含关系，那么选择区间小的那个环，舍弃区间大的那个环；
//		如果是相交且非包含关系，那么根据交集的大小进行抉择；
//		上面计算的目的是尽量将有交的环进行简化，减少切割点的个数（多环对应一个切割点）
	}//返回simplify_cycle_nodes 这是简化版的cycle_visiting_nodes
//	vector<vector<int>> cycle_samp_cutpos;
	int numcolors = colors.size()/dbg.size();
	string stdcolor(numcolors,'1');
	ssize_t numnodes = dbg.num_nodes();
	int kmer = dbg.k - 1;
	cutfile << "Cycle number:" << cycle_num << endl;
	for(int cy = 0;cy < cycle_num;cy++){
		vector<ssize_t> simplify_node = cycle_visiting_nodes[cy];
		cout << "cycle " << cy << "'s nodes:" << simplify_node << endl;
		//vector<ssize_t> simplify_node = simplify_cycle_nodes[cy];
		ssize_t cycle_startnode = simplify_node[0];
		ssize_t cycle_endnode = simplify_node[1];
		int startpos = nodepos[cycle_startnode];
		//int startpos = nodepos[cycle_startnode] - kmer + 1;
		//int endpos = nodepos[cycle_endnode];//can't set +kmer-1,eg.TATAT,the position of two repeat has overlap
		int endpos = nodepos[cycle_endnode]+kmer-1;
//		cout << "startnode " << cycle_startnode << ":" << startpos << endl;
//		cout << "endnode " << cycle_endnode << ":" << endpos << endl;
		vector<int> visited(numnodes,-1);
		map<int,int> pos2thick;
//		map<int,char> pos2bases;
		map<int,vector<char>> pos2bases;
		map<int,vector<string>> pos2color;//这个变量可以省略
//		map<int,string> pos2color;//这个变量可以省略
		map<int,int> pos2uniq;
		queue<ssize_t> nodelist;
		nodelist.push(cycle_startnode);
		while(!nodelist.empty()){
			ssize_t node = nodelist.front();
			int pos = nodepos[node];
			string label = dbg.node_label(node);
		//	cout << "node->pos" << node << "\t" << pos << "\t" << label << "\t" << getnodecolor(node) << endl;
			//if(!(pos2bases.find(pos))){
			vector<char> tmpbases = pos2bases[pos];
			if(tmpbases.empty()){
		//		cout << "1:node->pos " << node << "\t" << pos << "\t" << label << "\t" << getnodecolor(node) << endl;
				vector<char> tmpchar = {label[kmer - 1]};
				pos2bases[pos] = tmpchar;
			//	pos2bases[pos] = label[kmer - 1];
				pos2thick[pos] = 1;
				string tmpcolor1 = getnodecolor(node);
				tmpcolor1 = ColorCap({tmpcolor1,researched_samples});
				vector<string> tmpcolor = {tmpcolor1};
				pos2color[pos] = tmpcolor;;
//				pos2color[pos] = getnodecolor(node);
				pos2uniq[pos] = 1;	
			}
			else{
		//		cout << "2:node->pos " << node << "\t" << pos << "\t" << label << "\t" << getnodecolor(node) << "\t" << pos2bases[pos] << endl;
				vector<char> tmpbases = pos2bases[pos];
				int tmpnum = tmpbases.size();
				string tmpcolor = getnodecolor(node);
				tmpcolor = ColorCap({tmpcolor,researched_samples});
				vector<string> tmpcolors = pos2color[pos];
				int basesame = -1;
				for(int i = 0;i < tmpnum;i++){
					if(tmpbases[i] == label[kmer - 1]){
//						string newcolor = ColoCup({tmpcolor,tmpcolors[i]});
//						pos2color[pos] = newcolor;
						basesame = i;

						break;
					}
				}
				if(basesame == -1){
					tmpbases.push_back(label[kmer - 1]);
					pos2bases[pos] = tmpbases;
					tmpcolors.push_back(tmpcolor);
					pos2color[pos] = tmpcolors;
					pos2uniq[pos] = 0;
					pos2thick[pos]++;
				}else{
					string cupcolor = ColorCup({tmpcolors[basesame],tmpcolor});
					tmpcolors[basesame] = cupcolor;
					pos2color[pos] = tmpcolors;
				}
//				cout << "2:node->pos " << node << "\t" << pos << "\t" << label << "\t" << pos2color[pos] << "\t" << pos2bases[pos] << endl;
			}
			for(int x = 1;x < 5;x++){
				ssize_t nextnode = dbg.incoming(node,x);
				if(nextnode == -1||visited[nextnode] != -1 || nodepos[nextnode] > endpos){
					continue;
				}
				nodelist.push(nextnode);
				visited[nextnode] = 1;
			}
			for(int x = 1;x < 5;x++){
				ssize_t frontnode = dbg.outgoing(node,x);
				if(frontnode == -1 || visited[frontnode] != -1 || nodepos[frontnode] < startpos){
					continue;
				}
				nodelist.push(frontnode);
				visited[frontnode] = 1;
			}
			visited[node] = 2;
			nodelist.pop();
		}
		map<int,int> start2len;
		map<int,int> start2thick;
//		map<int,string> start2label;
		map<int,vector<string>> start2label;
		map<int,vector<string>> start2color;
//		int start = startpos;
//		int uniqlen = 0;
		string cutlabel = "";
		cutfile << "cycle:" << cy << endl;
		cutfile << "cout\ttype\tstar\tlength\tcutpos\tthick\trate\tlabel\tcolor" << endl; //type:1表示碱基优先；2表示thick优先
		int bestcutpos = -1;
		double bestvalue = 0.0;//rate = L/thick
		int cutid = 0;
//		for(int p = startpos;p <= endpos;p++){//在找切割点时，碱基一致优先
//		//	cout << "postion:" << p << "\tpos2uniq:" << pos2uniq[p] << "\tbase:" << pos2bases[p] << "\tpos2color:" << pos2color[p] << endl;
//			vector<string> tmpcolors = pos2color[p];
//			vector<char> tmpbases = pos2bases[p];
//			int tmpthick = tmpbases.size();
//			vector<string> cutlabels(tmpthick,"");
//			if(pos2uniq[p] == 1 && tmpcolors[0] == stdcolor){
//				if(uniqlen == 0){
//					start = p;
//					cutlabels[0] = "";
//				}			
//				uniqlen++;
//				cutlabels[0] = tmpbases[p] + cutlabels[0];
//				//cutlabel = pos2bases[p] + cutlabel;
//			}
//		//	else{
//			if(p == endpos || !(pos2uniq[p] == 1 && tmpcolors[0] == stdcolor)){
//				if(uniqlen == 0){continue;}
//				else{
////				if(pos2color[p] == stdcolor){
//					cutid++;
//					start2len[start] = uniqlen;
//					start2thick[start] = pos2thick[start];
//					start2label[start] = cutlabels;
//					double rate = start2len[start]/pos2thick[start];
//					cutfile << cutid << "\t" << 1 << "\t" << start << "\t" << start+uniqlen-1 << "\t" << uniqlen << "\t" << pos2thick[start] << "\t" << cutlabels[0] << "\t" << tmpcolors[0] << "\t" << rate << endl;
//					if(bestcutpos == -1){
//						bestcutpos = start;
//						bestvalue = uniqlen/pos2thick[start];
//					}
//					else{
//						if(pos2thick[start] == 1){
//							if(pos2thick[bestcutpos] == 1){
//								if(start2len[start] > start2len[bestcutpos]){
//									bestcutpos = start;
//									bestvalue = start2len[start];
//								}
//							}
//							else{
//								bestcutpos = start;
//								bestvalue = start2len[start];
//							}
//						}
//						else{
//							if(pos2thick[bestcutpos] == 1){continue;}
//							else{
//					//			double rate = start2len[start]/pos2thick[start];
//								if(rate > bestvalue){
//									bestcutpos = start;
//									bestvalue = rate;
//								}
//							}
//						}
//					}
//					uniqlen = 0;
//				}
//			}
//		}
//		int cutnum = start2len.size();
//		cutfile << "cut position number:" << cutnum << endl;
//		if(cutnum == 0){ //表示没有找到一个位点使得图在该位点上只有一个碱基（很大可能是没有对齐），接下来是碱基不一致时考虑thick优先				int newstart = 0;
			int newthick = pos2thick[startpos];
			int start = startpos;
			int uniqlen = 1;

			vector<string> cutcolors = pos2color[startpos];
			string poscolor = ColorCup(cutcolors);
			//vector<char> tmpbases = pos2bases[startpos];
			//int tmpthick = tmpbases.size();
			vector<string> cutlabels(4,"");
			//string cutlabel = "";
			for(int p1 = startpos;p1 <= endpos;p1++){
				int tmpthick = pos2thick[p1];
				vector<string> tmpcolors = pos2color[p1];
				string cupcolor = ColorCup(tmpcolors);
				cout << start << "\t" << p1 << "\t" << tmpthick << "\t" << newthick << "\t" << cupcolor << "\n";
					vector<char> tmpbases = pos2bases[p1];
				//if(tmpthick == newthick && cupcolor == stdcolor){
				if(tmpthick == newthick && poscolor == cupcolor){
					for(int i = 0;i < newthick;i++){
						string tmplabel = cutlabels[i];
						tmplabel = tmpbases[i] + tmplabel;
						cutlabels[i] = tmplabel;
						cout << tmplabel << "\t";
//						if(cutcolors[i] != tmpcolors[i]){
//							cout << "Meet error: same branch has different colors" << endl;
//						}
					}
					cout << endl;
					uniqlen++;					
				}
				if(!(tmpthick == newthick && poscolor == cupcolor) || p1 == endpos){
				//if(!(tmpthick == newthick && cupcolors == stdcolor) || p1 == endpos){
					cutid++;
//				else{					
					start2len[start] = uniqlen;
					start2thick[start] = newthick;
					//vector<string> 
					double a = uniqlen;
					double b = pos2thick[start];
					double rate = a/b;
//					double rate = uniqlen/pos2thick[start];
					for(int i = 0;i < newthick;i++){
						string finallabel = cutlabels[i];
						string finalcolor = cutcolors[i];
						int cutpos = (start + start + uniqlen - 1)/2;
				//		if(cupcolor == stdcolor)
						if(cupcolor == researched_samples){
							if(newthick == 1){
								cutfile << cutid << "\t" << 1 << "\t" << start << "\t"  << uniqlen << "\t" << cutpos << "\t" << newthick << "\t" << rate << "\t" << finallabel << "\t" << finalcolor << endl;
							}else{
								cutfile << cutid << "\t" << 2 << "\t" << start << "\t" << uniqlen << "\t" << cutpos << "\t" << newthick << "\t" << rate << "\t" << finallabel << "\t" << finalcolor  << endl;
							}
						}else{
							cutfile << cutid << "\t" << 3 << "\t" << start << "\t" <<  uniqlen << "\t" << cutpos << "\t" << newthick << "\t" << rate << "\t" << finallabel << "\t" << finalcolor << endl;
						}
						
					}
			//		cutfile << endl;
					if(bestcutpos == -1){
						bestcutpos = start;
						bestvalue = uniqlen/pos2thick[start];
					}
					else{
					//	double rate = uniqlen/pos2thick[start];
						if(bestvalue < rate){
							bestcutpos = start;
							bestvalue = rate;
						}
					}
					newthick = tmpthick;
					start = p1;
					uniqlen = 1;
					cutcolors = pos2color[start];
					poscolor = cupcolor;
					for(int i = 0;i < newthick;i++){
						cutlabels[i] = tmpbases[i];
					}
				}
			}	
//		}
		int final_cutpos = (bestcutpos+bestcutpos+start2len[bestcutpos]-1)/2;
		cout << "finally reference cut position:" << final_cutpos << "\treference rate:" << bestvalue<< endl;
		cutfile << "finally reference cut position:" << final_cutpos << "\treference rate:" << bestvalue<< endl;
//		vector<int> vecsamp_cutpos(numcolors,0);
//		for(ssize_t i = 0;i < numnodes;i++){
//			if(nodepos[i]){
//				if(nodepos[i] == final_cutpos){
//					string tmpcolor = getnodecolor(i);
//					cout << "find the cut node: " << i << "\tcolor: " << tmpcolor << endl;
//					for(int i1 = 0;i1 < numcolors;i1++){
//						if(tmpcolor[i1] == '1'){
//							map<ssize_t,int> sampnodepos = samp_nodepos[i1];
//							cout << "sample " << i1 << "'s nodepos: " << sampnodepos[i] << endl;
//							vecsamp_cutpos[i1] = sampnodepos[i];
//						}
//					}
//				}
//			}
//		}
//		cycle_samp_cutpos.push_back(vecsamp_cutpos);
	}
//	return cycle_samp_cutpos;
}

struct C_Su{
	vector<vector<ssize_t>> pairnodesets;
	vector<vector<int>> affiliationsets;
//	vector<int> nestlevelsets;
	vector<string> bubblecolors;
	map<ssize_t,int> nodepos;
	map<ssize_t,int> refnodepos;
	map<int,int> bub2gap;
	map<int,int> bub2var;
//	string researched_samples;
};

C_Su GetColoredSuperbubble(){
//C_Su GetColoredSuperbubble(int mystartnodepos){
	int t = getMilliCount();
//	myColCovPlus colcovplus = GetnewInfo();
//	vector<int> myplusnode = colcovplus.myplusnode;
//	vector<int> mycolors = colcovplus.mycolors;
//	vector<int> myplusnode = read_mystrand_from_file("strandfile.txt");
	int numcolors = colors.size()/dbg.size();
	cout << "The number of samples: " << numcolors <<endl;
	string stdcolor(numcolors,'1');
	string zerocolor(numcolors,'0');
	vector<ssize_t> noincomingnode;
	queue<ssize_t> supernodeset;
	C_Su colored_su;
	ssize_t numnodes = dbg.num_nodes();
	vector<int> myplusnode(numnodes,1);
	vector<vector<ssize_t>> pairnodesets;
	vector<vector<int>> affiliationsets;
	vector<string> bubblecolors;
//	vector<int> nestlevelsets(dbg.num_nodes(),1);
	int buborder = 0;
	map<int,bool> bubend;
	map<ssize_t,int> nodepos;
	map<ssize_t,int> refnodepos;
//	ofstream nestlevelfile("mynestlevel.txt");
	ofstream pairsfile("mypairnodes.txt");
//	ofstream indelfile("myindelinfo.txt");
	ofstream posfile("mynodeposinfo.txt");
	ofstream msafile("msaallsamples.txt");
	ofstream varfile("myvar.txt");
//	ofstream cutposfile("mycutposinfo.txt");
//	ofstream affiliationfile("myaffiliation.txt");

//	string zerolabel = dbg.node_label(0);
//	int kmer = zerolabel.size();
	int kmer = dbg.k-1;
	map<ssize_t,int> node2gap;
	map<ssize_t,string> node2gapcolor;
	int endstringlen = dbgendstring.size();
	string endlabel = dbgendstring.substr(endstringlen - kmer,kmer);
	cout << "endnodelabel:" << endlabel << endl;
	string startlabel = dbgstartstring.substr(1,kmer);
	ssize_t startnode = -1;
	string refstring;
	for(ssize_t node = 0;node < numnodes;node++){
//		if(dbg.outdegree(node) == 0){
//			cout << node << endl;
//		string tmpnodelabel = dbg.node_label(node);
//			cout << tmpnodelabel << endl;
//			if(getnodecolor(node) != stdcolor){
//				cout << "start node does not contain all samples\n";
//			}
//		}
//		else{continue;}
		string tmpnodelabel = dbg.node_label(node);
		if(tmpnodelabel == endlabel){
			noincomingnode.push_back(node);
			supernodeset.push(node);
			startnode = node;
			cout << "find start node: " << node << "\tpos:" << nodepos[node] << "\tseq:" << dbg.node_label(node) << endl;
			if(getnodecolor(node) != stdcolor){
				cout << "start node does not contain all samples\n";
				break;
			}
			nodepos[node] = 1;
			refnodepos[node] = 1;
			refstring = dbg.node_label(node);
//			posfile << node << "\t" << nodepos[node] << endl;
			break;
		}
	}

//	int noincomingnodenum = noincomingnode.size();
//	cout << "There are " << noincomingnodenum << " no incoming nodes" << endl;
//	cout << "startnode:" << startnode << endl;
	if(supernodeset.empty()){
		cout << "Can not find start node" << endl;
		exit(0);
	}
	vector<int> visited(dbg.num_nodes(),0);
	vector<int> sink_visited(dbg.num_nodes(),0);
	vector<int> visited_order(dbg.num_nodes(),0);
	vector<ssize_t> sourcenodeset;
	vector<vector<string>> sourcenodecolors;
	vector<string> bubcolors;
	vector<vector<ssize_t>> bubpairs;
	
//	supernodeset.push(startnode);
	bool getend = false;
	vector<vector<ssize_t>>	cycle_visiting_nodes;
	vector<ssize_t> cycle_visiting_node;
//	vector<map<ssize_t,int>> samp_nodepos(numnodes,nodepos);//
	map<ssize_t,int> visiting_nodes;
	visiting_nodes[startnode] = 1;
	int cycle_count = 0;
	map<ssize_t,int> cycle_startnodes;
	int maxnodepos = 0;
	int maxrefpos = 0;
	while(getend == false){
		visiting_nodes.erase(visiting_nodes.begin(),visiting_nodes.end());
		while(!supernodeset.empty()){
			ssize_t node = supernodeset.front();
	//		ssize_t fathernode = node;
			cout << "\nsupernode:" << node << endl;
			cout << "supernodelist's length: " << supernodeset.size() << endl;
				if(visiting_nodes[node]){
					visiting_nodes.erase(node);
						cout << "test if erase 0? " << node << "\t" << visiting_nodes[node] << endl;
				}
			string nodecolor = getnodecolor(node);
			nodecolor = ColorCap({nodecolor,researched_samples});
			if(researched_samples == nodecolor){
				cout << "stdcolor node :" << node << endl;
			}
	//		if(node == 25303)
	//			testbool = false;
	//		if(testbool == false){
	//			if(getnodecolor(node) == stdcolor){
	//				cout << "stdcolor node:" << node << endl;
	//			}
	//		}
			int intimes = 0;
			for(unsigned long x = 1;x < dbg.sigma+1;x++){
				ssize_t nextnode = dbg.incoming(node,x);
				if(nextnode == -1 || myplusnode[nextnode] == 0)
					continue;
				
				if(cycle_startnodes[nextnode] && cycle_startnodes[nextnode] > 0){
					int tmpcycleid = cycle_startnodes[nextnode];
					vector<ssize_t> tmpcycle_visiting_node = cycle_visiting_nodes[tmpcycleid-1];
							cout << node << "<-" << nextnode << endl;
					cout << "add cycle endnode:" << tmpcycle_visiting_node << "<-" << node << endl;
					tmpcycle_visiting_node.push_back(node);
					cycle_visiting_nodes[tmpcycleid-1] = tmpcycle_visiting_node;
				}
				if(visited[nextnode] == 2){
					cout << "meet a node whose station is 2:" << nextnode << endl;
					if(cycle_visiting_node.empty()){
						cout << "error" << endl;
					}else{
				//		cycle_visiting_node.push_back(nextnode);
					}
					continue;
	//				cout << node << "->" << nextnode << " visit ordering has mistakes" << endl;
				}
				string nextnodecolor = getnodecolor(nextnode);
				string actualcolor = ColorCap({nextnodecolor,researched_samples});
				if(actualcolor == zerocolor)
					continue;
				ssize_t nextedge = dbg.incoming_edge(node,x);
				string nextedgecolor = getedgecolor(nextedge);
				nextedgecolor = ColorCap({nextedgecolor,researched_samples});
				if(nextedgecolor[numcolors-1] == '1'){
				//if(actualcolor[numcolors-1] == '1'&& nodecolor[numcolors-1] == '1')
					refnodepos[nextnode] = refnodepos[node] + 1;
					refstring = base[x] + refstring;
				}
	//			if(nodepos[nextnode] && abs(nodepos[node]-nodepos[nextnode]) > 1){
	//			if(nodepos[nextnode]){
	//				if(nodepos[node]-nodepos[nextnode] > -1){
	//					nodepos[nextnode] = nodepos[node] + 1;
	//					for(int i = 0;i < numcolors;i++){
	//						if(nextedgecolor[i] == '1'){
	//							map<ssize_t,int> tmpnodepos = samp_nodepos[i];
	//							tmpnodepos[nextnode] = nodepos[node] + 1;
	//							samp_nodepos[i] = tmpnodepos;
	//						}
	//					}
	//				}
	//				indelfile << node << ":" << nextnode << "\t" << nodepos[node] << "\t" << nodepos[nextnode] << endl;
	//			}
//					for(int i = 0;i < numcolors;i++){
//						if(nextedgecolor[i] == '1'){
//							map<ssize_t,int> tmpnodepos = samp_nodepos[i];
//							tmpnodepos[nextnode] = tmpnodepos[node] + 1;
//							samp_nodepos[i] = tmpnodepos;
//						}
//					}
				if(!nodepos[nextnode]){
					nodepos[nextnode] = nodepos[node]+1;
					node2gapcolor[nextnode] = nextedgecolor;
				}
				if(nodepos[nextnode] && nodepos[nextnode] != nodepos[node] + 1){
					if(!node2gap[nextnode]){
						if(nodepos[nextnode] < nodepos[node] + 1){
							node2gap[nextnode] = nodepos[node] - (nodepos[nextnode] - 1);
						}else if(nodepos[nextnode] > nodepos[node] + 1){
							node2gap[nextnode] = nodepos[nextnode] -1 - nodepos[node];
							node2gapcolor[nextnode] = nextedgecolor;
						}
					}else{
						int tmpmax = -1;
						int tmpmin = numnodes;
						string tmpcolor = zerocolor;
						for(int i = 1;i < 5;i++){
							ssize_t tmpoutnode = dbg.outgoing(nextnode,i);
							if(tmpoutnode == -1)
								continue;
			//				string tmpoutcolor = ColorCap({getnodecolor(tmpoutnode),researched_samples});
							if(!nodepos[tmpoutnode]){
								continue;
							}
							for(int j = 1;j < 5;j++){
								ssize_t node1 = dbg.incoming(tmpoutnode,j);
								if(node1 == nextnode){
									ssize_t edge1 = dbg.incoming_edge(tmpoutnode,j);
									tmpcolor = ColorCap({getedgecolor(edge1),researched_samples});
									break;
								}
							}
							if(tmpcolor == zerocolor){
								continue;
							}
							
							if(tmpmax < nodepos[tmpoutnode])
								tmpmax = nodepos[tmpoutnode];
							if(tmpmin == nodepos[tmpoutnode] && tmpcolor[numcolors - 1] == '1')
								node2gapcolor[nextnode] = tmpcolor;
							if(tmpmin > nodepos[tmpoutnode]){
								tmpmin = nodepos[tmpoutnode];
								node2gapcolor[nextnode] = tmpcolor;
							}
						}
						node2gap[nextnode] = tmpmax-tmpmin;
					}
					if(nodepos[nextnode] < nodepos[node] + 1){
						nodepos[nextnode] = nodepos[node]+1;
					}
				}
//				posfile << nextnode  << "\t" << nodepos[nextnode] << endl;
			//	if(nextnode == 33745){
			//		cout << nextnode  << "\t" << nodepos[nextnode] << "\t" << node << "\t" << nodepos[node] << endl;
			//	}
			
	//			if(nodepos[node].empty()){nodepos[nextnode] = nodepos[node]+1;}
	//			else{
	//				if(nodepos[nextnode] < nodepos[node] + 1){
	//					nodepos[nextnode] = nodepos[node]+1;
	//				}
	//			}
	//			bool bridge = false;
	//			if(nextnodecolor == stdcolor){
	//				nestlevelsets[nextnode] = 1;
	//				bridge = true;
	//			}
	//			bool getstart = false;
			//	while((dbg.outdegree(nextnode) == 1 && dbg.indegree(nextnode) == 1)||getstart == false){
	//			bool nochild = true;
				while(dbg.outdegree(nextnode) == 1 && dbg.indegree(nextnode) == 1){
	//				cout << "->" << nextnode;
	//				if(bridge == true ||nestlevelsets[nextnode] == 0)
	//					nestlevelsets[nextnode] = 1;
					visited[nextnode] = 2;
					int incomingid = -1;
					for(unsigned long x2 = 1;x2 < dbg.sigma+1;x2++){
						ssize_t nextnode2 = dbg.incoming(nextnode,x2);
						if(nextnode2 == -1)
							continue;
						if(cycle_startnodes[nextnode2] && cycle_startnodes[nextnode2] > 0){
							int tmpcycleid = cycle_startnodes[nextnode2];
							vector<ssize_t> tmpcycle_visiting_node = cycle_visiting_nodes[tmpcycleid-1];
							//cout << nextnode << "<-" << nextnode2 << endl;
							cout << "add cycle endnode:" << tmpcycle_visiting_node << "<-" << nextnode << endl;
							tmpcycle_visiting_node.push_back(nextnode);
							cycle_visiting_nodes[tmpcycleid-1] = tmpcycle_visiting_node;
						}
						if(visited[nextnode2] == 2)
							continue;
					//	ssize_t nextedge2 = dbg.incoming_edge(nextnode,x2);
					//	string nextedgecolor2 = getedgecolor(nextedge2);
						//nextedgecolor2 = ColorCap({nextedgecolor2,researched_samples});
						if(nextedgecolor[numcolors-1] == '1'){
						//if(actualcolor[numcolors-1] == '1')
							refnodepos[nextnode2] = refnodepos[nextnode] + 1;
							refstring = base[x2] + refstring;
						}
				//		string actualcolor2 = ColorCap({nextnodecolor,researched_samples});
				//		if(actualcolor2 == zerocolor)
				//			continue;
	//					if(nodepos[nextnode2]){
	//						if(nodepos[nextnode]-nodepos[nextnode2] > -1){
	//							nodepos[nextnode2] = nodepos[nextnode] + 1;
	//						}
	//						indelfile << nextnode << "\t" << nextnode2 << "\t" << nodepos[nextnode] << "\t" << nodepos[nextnode2] << endl;
	//					}
//							for(int i = 0;i < numcolors;i++){
//								if(nextedgecolor[i] == '1'){
//									map<ssize_t,int> tmpnodepos = samp_nodepos[i];
//									tmpnodepos[nextnode2] = tmpnodepos[nextnode] + 1;
//									samp_nodepos[i] = tmpnodepos;
//								}	
//							}
						if(!nodepos[nextnode2]){
							nodepos[nextnode2] = nodepos[nextnode]+1;
							node2gapcolor[nextnode2] = nextedgecolor;
						}
						if(nodepos[nextnode2] && nodepos[nextnode2] != nodepos[nextnode] + 1){
							if(!node2gap[nextnode2]){
								if(nodepos[nextnode2] < nodepos[nextnode] + 1){
									node2gap[nextnode2] = nodepos[nextnode] - (nodepos[nextnode2] - 1);
								}else if(nodepos[nextnode2] > nodepos[nextnode] + 1){
									node2gap[nextnode2] = nodepos[nextnode2] -1 - nodepos[nextnode];
									node2gapcolor[nextnode2] = nextedgecolor;
								}
							}else{
								int tmpmax = -1;
								int tmpmin = numnodes;
								string tmpcolor = zerocolor;
								for(int i = 1;i < 5;i++){
									ssize_t tmpoutnode = dbg.outgoing(nextnode2,i);
									if(tmpoutnode == -1)
										continue;
									if(!nodepos[tmpoutnode]){
										continue;
									}
					//				string tmpoutcolor = ColorCap({getnodecolor(tmpoutnode),researched_samples});
									for(int j = 1;j < 5;j++){
										ssize_t node1 = dbg.incoming(tmpoutnode,j);
										if(node1 == nextnode2){
											ssize_t edge1 = dbg.incoming_edge(tmpoutnode,j);
											tmpcolor = ColorCap({getedgecolor(edge1),researched_samples});
											break;
										}
									}
									if(tmpcolor == zerocolor){
										continue;
									}
									
									if(tmpmax < nodepos[tmpoutnode])
										tmpmax = nodepos[tmpoutnode];
									if(tmpmin == nodepos[tmpoutnode] && tmpcolor[numcolors - 1] == '1')
										node2gapcolor[nextnode2] = tmpcolor;
									if(tmpmin > nodepos[tmpoutnode]){
										tmpmin = nodepos[tmpoutnode];
										node2gapcolor[nextnode2] = tmpcolor;
									}
								}
								node2gap[nextnode2] = tmpmax-tmpmin;
							}
							if(nodepos[nextnode2] < nodepos[nextnode] + 1){
								nodepos[nextnode2] = nodepos[nextnode]+1;
							}
						}
							//indelfile << nextnode << "\t" << nextnode2 << "\t" << nodepos[nextnode] << "\t" << nodepos[nextnode2] << endl;
							
				//		}
//						posfile << nextnode2 << "\t" << nodepos[nextnode2] << endl;
//				if(nextnode2 == 33745){
//					cout << nextnode2  << "\t" << nodepos[nextnode2] << "\t" << nextnode << "\t" << nodepos[nextnode] << endl;
//				}
						if(dbg.indegree(nextnode2) == 0){
							nextnode = nextnode2;
							break;
						}
//						if(myplusnode[nextnode2] == 0){
//							visited[nextnode2] = 2;
//							continue;
//						}
	//					fathernode = nextnode;
						nextnode = nextnode2;
	//					nochild = false;
						incomingid++;
						break;
					}
					if(incomingid == -1){break;}
	//				cout << nextnode << ":" << dbg.outdegree(nextnode) << ":" << dbg.indegree(nextnode) << endl;
				}
	//		if(nochild == true){continue;}
	//			cout << "->" << nextnode << endl;
	//			cout << "indegree:" << dbg.indegree(nextnode) << " outdegree:" << dbg.outdegree(nextnode) << " visited:" << visited[nextnode] << endl;
	//			if(nextnode == 86567){cout << "get the end of DBG";}
				bool innodevisited = true;
	//			if(getnodecolor(nextnode) == stdcolor){
	//				cout << "stdcolor node2 :" << nextnode << endl;
	//			}
				if(dbg.outdegree(nextnode) > 1){//sink node
					string circlecol = zerocolor;
					if(sink_visited[nextnode] == 1){
						continue;
					}
					for(unsigned long x2 = 1;x2 < dbg.sigma+1;x2++){
						ssize_t frontnode = dbg.outgoing(nextnode,x2);
						if(frontnode == -1 || myplusnode[frontnode] == 0)
							continue;
						string tmpcolor = getnodecolor(frontnode);
						string actualcolor2 = ColorCap({tmpcolor,researched_samples});
						if(actualcolor2 == zerocolor)
							continue;
	
	//					string tmpnodecolor = getnodecolor(frontnode);
	//					int colcmp = StringCMP(circlecol,tmpnodecolor);
	//					if(colcmp == 1 ||(colcmp == 2 && tmpnodecolor != zerocolor)){
	//						cerr << "meet circle at " << nextnode << endl; 
	///					}
	//					else{
	//						circlecol = ColorCup({circlecol,tmpnodecolor});
	//					}
	
						if(visited[frontnode] == 0){
							innodevisited = false;
							visiting_nodes[nextnode] = 1;
							cout << "visiting node: " << nextnode << endl;
							break;
						}
					}
					if(innodevisited == true){//means all the incoming nodes of sink node have been visited
						if(visited[nextnode] < 1){//0 or -1
							supernodeset.push(nextnode);
							cout << "push a sink node " << nextnode << endl;
							visited[nextnode] = 1;
						}
						if(visiting_nodes[nextnode] && visited[nextnode] > 0){
							visiting_nodes.erase(nextnode);
							cout << "test if erase?" << nextnode << "\t" << visiting_nodes[nextnode] << endl;
						}
						string nextnodecolor = getnodecolor(nextnode);
						string actualcolor2 = ColorCap({nextnodecolor,researched_samples});
						if(actualcolor2 == zerocolor)
							continue;
	//					cout << "sink supernode:" << nextnode << "\tcolor:" << nextnodecolor << endl;
	//					if(nextnodecolor == zerocolor)
	//						continue;
	//					if(nextnodecolor != stdcolor)
	//						continue;
	///////////// To find the matched source node ////
	//		if(testbool == false){
	//			if(nextnodecolor == stdcolor){
	//				cout << "stdcolor node:" << nextnode << endl;
	//			}
	//		}
						vector<string> inedgecolor;//store the incoming edges' color information of sink node
						int incount = 0;
						for(unsigned long x3 = 1;x3 < dbg.sigma+1;x3++){
							ssize_t frontnode = dbg.outgoing(nextnode,x3);
							if(frontnode == -1 || myplusnode[frontnode] == 0)
								continue;
							for(unsigned long x4 = 1;x4 < dbg.sigma+1;x4++){
								ssize_t tmpnextnode = dbg.incoming(frontnode,x4);
								if(tmpnextnode == -1 || tmpnextnode != nextnode)
									continue;
								ssize_t tmpedge = dbg.incoming_edge(frontnode,x4);
								string edgecolor = getedgecolor(tmpedge);
								edgecolor = ColorCap({edgecolor,researched_samples});
								if(edgecolor == zerocolor)
									continue;
								//if(edgecolor == zerocolor)
								//	continue;
					//				edgecolor = getnodecolor(nextnode);
								if(inedgecolor.empty()){
	//								cout << "push back inedgecolor:" << edgecolor << endl;
									//inedgecolor.push_back(edgecolor);
									inedgecolor.push_back(edgecolor);
									incount++;
								}
								else{
	//								cout << "push back inedgecolor2:" << edgecolor << endl;
									string cupstr = ColorCup(inedgecolor);
							//		string strcol2 = getedgecolor(tmpedge);
							//		string strcol1 = inedgecolor[i];
									int strcmp = StringCMP(cupstr,edgecolor);
	//								cout << "sink strcmp:" << strcmp << " str1: " << cupstr << "\tstr2: " << edgecolor << endl;
									if(strcmp == 0){
										inedgecolor.push_back(edgecolor);
										incount++;
									}
									else{
										string nodelabel = dbg.node_label(nextnode);
	//									cout << nextnode << ":" << nodelabel << endl;
										if(nodelabel[0] != '$')
											cout << nextnode << " has back edge" << endl;
										
	
	
									}
								}
								break;
							}
						}
						if(incount > 1){
							cout << "find a sink node " << nextnode << ";start to match its source node"  << endl;
	//##########	########  start to match source node  ####################
							int sourcenum = sourcenodeset.size();
	//						int sourcecolornum = sourcenodecolors.size();
	//						cout << sourcecolornum << endl;
							cout << "There are " << sourcenum << " source nodes in the source node set" << endl;
							if(sourcenum == 0){
								cout << "sink node: " << nextnode << " can't find the matched source node" << endl;
	//							cout << "node color: " << getnodecolor(nextnode) << endl;
								continue;
							}
							int insize = inedgecolor.size();
							string incolorcup = ColorCup(inedgecolor);
							bool sinkjudge = false;
							for(int i = sourcenum - 1;i >= 0;i--){
								vector<string> outedgecolor = sourcenodecolors[i];
	//							bool colorsetcmp = true;
								int outsize = outedgecolor.size();
	//							cout << "insize:" << insize << "\toutsize:" << outsize << endl;
								if(outsize < 2 || insize < 2){cout << "source node or sink node colorset meets error!" << endl;}
	//							cout << inedgecolor[0] << ":" << inedgecolor[1] << endl;
	//							cout << outedgecolor[0] << ":" << outedgecolor[1] << endl;
								string outcolorcup = ColorCup(outedgecolor);
	//							cout << "incolorcup:" << incolorcup << "\toutcolorcup:" << outcolorcup << endl;
								int incapcount = 0;
								int outcapcount = 0;
								for(int i1 = 0;i1 < insize;i1++){
									string inedge = inedgecolor[i1];
									vector<string> tmpcolor = {inedge,outcolorcup};
									string capstr = ColorCap(tmpcolor);
									if(capstr != zerocolor)
										incapcount++;
								}
								for(int i2 = 0;i2 < outsize;i2++){
									string outedge = outedgecolor[i2];
									vector<string> tmpcolor = {outedge,incolorcup};
									string capstr = ColorCap(tmpcolor);
									if(capstr != zerocolor)
										outcapcount++;
								}
								ssize_t sourcenode = sourcenodeset[i];
								if(sourcenode == nextnode){
									cout << "same source and sink node: " << nextnode << endl;
				//					sourcenodeset.erase(sourcenodeset.begin()+i);
				//					sourcenodecolors.erase(sourcenodecolors.begin()+i);
				//					continue;
								}
	//							cout << "incapcount:" << incapcount << "\toutcapcount:" << outcapcount << endl;
								if(incapcount > 1 && outcapcount > 1){
									//string sourcecolor = getnodecolor(sourcenode);
									//string sinkcolor = getnodecolor(nextnode);
									string sourcecolor = outcolorcup;
									string sinkcolor = incolorcup;
									vector<string> allcolorstr = {incolorcup,outcolorcup};
	//								cout << "sourcecolor:" << sourcecolor << " \tsinkcolor:" << sinkcolor << endl;
									string commoncolor = ColorCap(allcolorstr);
	//								cout << "commoncolor:" << commoncolor << endl;
									if(sourcecolor == sinkcolor|| commoncolor == sourcecolor){// need to delete sourcenode
	//									cout << "find matched nodes(delete):" << sourcenode << ":" << nextnode << "\t" << commoncolor << endl;
										
										sourcenodeset.erase(sourcenodeset.begin()+i);
										sourcenodecolors.erase(sourcenodecolors.begin()+i);
										int tmpsourcenum = sourcenodeset.size();
										if(tmpsourcenum == 0)
											bubend[buborder] = true;
									}
									if(sourcecolor == sinkcolor || commoncolor == sinkcolor){
										sinkjudge = true;
									}
								//	else{
								//		cout << "find matched nodes(undeleted):" << sourcenode << ":" << nextnode << "\t" << commoncolor << endl;
	
								//	}
									buborder++;
	//								pairsfile << buborder << "\t" << sourcenode << ":" << nextnode << "\t" << sourcecolor << ":" << sinkcolor << "\t" << commoncolor << endl;
									//vector<ssize_t> pairs = {sourcenode,nextnode};
									vector<ssize_t> pairs = {nextnode,sourcenode};
									bubpairs.push_back(pairs);
									bubcolors.push_back(commoncolor);
									cout << "find a cSupB:" << nextnode << ":" << sourcenode << endl;
								}
								if(sinkjudge == true){
									break;
								}
							}
							sink_visited[nextnode] = 1;
							//if(incolorcup == stdcolor)
							if(incolorcup == researched_samples){
								int sournum = sourcenodeset.size();
								if(sournum > 0){cout << "delete source node has problem" << endl;}
							}
						}
						
					}
					else{
						continue;
					}
				}
				if(dbg.indegree(nextnode) > 1 && innodevisited == true){//source node
					//bool issource = true;
					int outcount = 0;
		//			if(dbg.outdegree(nextnode) < 1){
					if(visited[nextnode] < 1){
						supernodeset.push(nextnode);
						visited[nextnode] = 1;
					}
					cout << "push a source node " << nextnode << endl;
					if(visiting_nodes[nextnode] && visited[nextnode] > 0){
						visiting_nodes.erase(nextnode);
							cout << "test if erase 2?" << nextnode << "\t" << visiting_nodes[nextnode] << endl;
					}
				//	visiting_nodes[nextnode] = 0;
					string nextnodecolor = getnodecolor(nextnode);
	//				cout << "source supernode:" << nextnode << "\tcolor:" << nextnodecolor <<  endl;
					string actualcolor2 = ColorCap({nextnodecolor,researched_samples});
					if(actualcolor2 == zerocolor)
						continue;
		//			if(nextnodecolor == zerocolor)
		//				continue;
	//				if(nextnodecolor != stdcolor)
	//					continue;
	//		if(testbool == false){
	//			if(nextnodecolor == stdcolor){
	//				cout << "stdcolor node:" << nextnode << endl;
	//			}
	//		}
					vector<string> outedgecolor;
					for(unsigned long x2 = 1;x2 < dbg.sigma+1;x2++){ // need to test if nextnode is a source node based on colors
						ssize_t nextnode2 = dbg.incoming(nextnode,x2);
						if(nextnode2 == -1)
							continue;
						for(unsigned long x3 = 1;x3 < dbg.sigma+1;x3++){
							ssize_t frontnode = dbg.outgoing(nextnode2,x3);
							if(frontnode == -1 || frontnode != nextnode)
								continue;
							ssize_t tmpedge = dbg.incoming_edge(nextnode,x2);
							string edgecolor = getedgecolor(tmpedge);
							edgecolor = ColorCap({edgecolor,researched_samples});
							if(edgecolor == zerocolor)
								continue;
						//		edgecolor = getnodecolor(nextnode);
							if(outedgecolor.empty()){
	//							cout << "push back outedgecolor:" << edgecolor << endl;
								outedgecolor.push_back(edgecolor);
								outcount++;
							}
							else{
								string cupstr = ColorCup(outedgecolor);
							//	string strcol2 = getedgecolor(tmpedge);
								int strcmp = StringCMP(cupstr,edgecolor);
	//							cout << "source strcmp:" << strcmp << " str1:" << cupstr << "\tstr2:" << edgecolor << endl;
								if(strcmp == 0){
									outedgecolor.push_back(edgecolor);
	//								cout << "push back outedgecolor2:" << edgecolor << endl;
									outcount++;
								}	
								else{
									string nodelabel = dbg.node_label(nextnode);
									cout << nextnode << ":" << nodelabel << endl;
									if(nodelabel[0] != '$')
										cout << nextnode << " has back edge" << endl;
	
	
								}
							}
							break;
						}
					}
					if(outcount > 1){ // nextnode is a source node
						sourcenodeset.push_back(nextnode);
						sourcenodecolors.push_back(outedgecolor);
					}
				}
//				if(dbg.indegree(nextnode) == 0 || getnodecolor(nextnode) == zerocolor){//end node
//					getend = true;
//					visited[nextnode] = 2;
//					cout << "Getting end,but this is not the result we want to see.You need to modify(change or extand) the head of the sequences" << endl;
//					maxnodepos = nodepos[nextnode];
//					continue;
//				}
				if(dbg.indegree(nextnode) == 0 || getnodecolor(nextnode) == zerocolor || dbg.node_label(nextnode) == startlabel){//end node
					getend = true;
					visited[nextnode] = 2;
					maxnodepos = nodepos[nextnode];
					maxrefpos = refnodepos[nextnode];
					continue;
				}
	//			if(visited[nextnode] == 0){
	//				visiting_nodes[nextnode] = 1;
	//			}
				if(visited[nextnode] > 0){
					if(visiting_nodes[nextnode]){
						visiting_nodes.erase(nextnode);
							cout << "test if erase 3?" << nextnode << "\t" << visiting_nodes[nextnode] << endl;
					}
				//	visiting_nodes[nextnode] = 0;
				}
				intimes++;
			}
			visited[node] = 2;
			supernodeset.pop();
		}
		
	//	if(!cycle_visiting_node.empty()){
	//		cycle_visiting_nodes.push_back(cycle_visiting_node);
	//	}
		cout << "getend: " << getend << endl;
		if(getend == false){
			cycle_count++;
			map<ssize_t,int>::iterator iter;
			vector<string> colorset;
			vector<ssize_t> actual_visiting_nodes;
			for(iter = visiting_nodes.begin();iter != visiting_nodes.end();iter++){
				ssize_t tmpnode = iter->first;
				int state = iter->second;
				if(state == 1){
					actual_visiting_nodes.push_back(tmpnode);
					string tmpcolor = getnodecolor(tmpnode);
					tmpcolor = ColorCap({tmpcolor,researched_samples});
					colorset.push_back(tmpcolor);
					cout << tmpnode << ":" << tmpcolor << endl;
				//int = ier
				}
			}
			int visitinglen = actual_visiting_nodes.size();//这里len>=2，否则存在问题
			cout << "There are " << visitinglen << " visiting nodes:" << actual_visiting_nodes << endl;
	//		for(int i = 0;i < visitinglen;i++){
	//			ssize_t tmpnode = visiting_nodes[i];
	//			cout << tmpnode << endl;
	//			string tmpcolor = getnodecolor(tmpnode);
	//			colorset.push_back(tmpcolor);
	//		}
	//
			string common_color = ColorCap(colorset);
			cout << "common color: " << common_color << endl;
			if(common_color == zerocolor){
				if(visitinglen == 2){
					cout << "common==zerocolor && visiting_length == 2 meets error" << endl;
				}
				else{
					//说明不止有一个环，这里暂时选择位置最靠前的node,这时可能会导致一些错误情况
					for(int i = 1;i < visitinglen;i++){
						if(nodepos[actual_visiting_nodes[0]]>nodepos[actual_visiting_nodes[i]]){
							ssize_t tmpnode = actual_visiting_nodes[0];
							actual_visiting_nodes[0] = actual_visiting_nodes[i];
							actual_visiting_nodes[i] = tmpnode; 
						}
					}
					ssize_t tmpstartnode = actual_visiting_nodes[0];
					//visited[tmpstartnode] = 1;
					cycle_visiting_node.clear();
					cycle_visiting_node.push_back(tmpstartnode);
					cout << "tmpstartnode:" << tmpstartnode << endl;
	//				cycle_visiting_nodes.push_back(cycle_visiting_node);
					supernodeset.push(tmpstartnode);
					visited[tmpstartnode] = 1;
				}
			}
			else{
	//			ssize_t startnode2 =actual_visiting_nodes[0];//表示最靠前的node
				ssize_t pos = actual_visiting_nodes[0];
	//			string commoncolor = ColorCap({getnodecolor(pos),getnodecolor(startnode2)});
				int select_samp = 0;
				for(int i1 = 0;i1 < numcolors;i1++){
					if(common_color[i1] == '1'){
						select_samp = i1;
						break;
					}
				}
				cout << "selected sample id: " << select_samp << endl;
				bool endsearch = false;
				ssize_t tmpstartnode = pos;
				while(endsearch == false){
					bool findfront = false;
					for(int x = 1;x < 5;x++){
						ssize_t frontnode = dbg.incoming(pos,x);
						if(frontnode == -1){
							continue;
						}
						string tmpcolor = getnodecolor(frontnode);
						tmpcolor = ColorCap({tmpcolor,researched_samples});
	//					cout << "tmpcolor: " << tmpcolor << endl;
						if(tmpcolor[select_samp] == '0'){continue;}
						pos = frontnode;
						findfront = true;
						if(std::find(actual_visiting_nodes.begin(),actual_visiting_nodes.end(),pos) != actual_visiting_nodes.end()){
					//	if(actual_visiting_nodes[pos] == 1){
							tmpstartnode = pos;
						}
						break;
					}
					if(findfront == false){
						endsearch = true;
					}
				}
				cout << "find the tmpstartnode: " << tmpstartnode << endl;
				cycle_startnodes[tmpstartnode] = cycle_count;
	//			visited[tmpstartnode] = 1;
				cycle_visiting_node.clear();
				cycle_visiting_node.push_back(tmpstartnode);
				cycle_visiting_nodes.push_back(cycle_visiting_node);
	//			cycle_visiting_nodes.push_back(cycle_visiting_node);
				supernodeset.push(tmpstartnode);
				visited[tmpstartnode] = 1;
			}		
		}
		visiting_nodes.erase(visiting_nodes.begin(),visiting_nodes.end());
	}
//	for(int i = 0;i < numcolors;i++){
//		map<ssize_t,int> tmpnodepos = samp_nodepos[i];
//		cout << "sample:" << i << "\tpostion:" << tmpnodepos[12509] << endl;;
//	}
	int cycle_num = cycle_visiting_nodes.size();		
	if(cycle_num == 0){
		cout << "No cycle in the cdbg" << endl;
	}
	else{
		cout << "find the cycle" << endl;
		//vector<vector<int>> cycle_samp_cutpos = Get_Cut_Pos(nodepos,samp_nodepos,cycle_visiting_nodes);
		Get_Cut_Pos(nodepos,cycle_visiting_nodes);
//		int cutnum = cycle_samp_cutpos.size();
//		for(int i0 = 0;i0 < cutnum;i0++){
//			vector<int> samp_cutpos = cycle_samp_cutpos[i0];
//			for(int i1 = 0;i1 < numcolors - 1;i1++){
//				cutposfile << samp_cutpos[i1] << ":";
//			}
//			cutposfile << samp_cutpos[numcolors - 1] << endl;
//		}
		cout << "To cut the the sequences in the Linux environment" << endl;
		exit(0);
	}
	cout << "Finished pairnodes searching work" << endl;
	ofstream reffile("myref.txt");
	reffile << refstring << endl;
	if(maxnodepos <= 0){
		cout << "max nodepos is error" << endl;
		exit(0);
	}
	map<ssize_t,int>::iterator it;
	cout << "max node position: " << maxnodepos << endl;
	for(it = nodepos.begin();it != nodepos.end();it++){
		ssize_t mynode = it->first;
		int mynodepos = it->second;
		if(refnodepos[mynode]){
			//int refnewpos = REF_LEN - refnodepos[mynode] - kmer + 2;
			int refnewpos = maxrefpos - refnodepos[mynode] + 1;
			refnodepos[mynode] = refnewpos;
		}else{
			refnodepos[mynode] = -1;
		}
		int newpos = maxnodepos + 1 - mynodepos;
		if(newpos < 0){
			cout << "max nodepos is error2" << endl;
			exit(0);

		}
		nodepos[mynode] = newpos;
		posfile << mynode << "\t" << newpos << "\t" << refnodepos[mynode] << "\t" << node2gap[mynode] << endl;
	}
//	for(ssize_t i = 0;i < dbg.num_nodes();i++){
//		posfile << i << "\t" << nodepos[i] << endl;
//	}
	int bubnum = bubcolors.size();
	cout << "There " << bubnum << " bubble found" << endl;
	map<int,int> bubfather;//store the father bubble of the bubble
	map<int,vector<int>> bubchild;//store the children bubble of the bubble
	for(int bub = 0;bub < bubnum;bub++){
		string color1 = bubcolors[bub];
//		cout << bub+1 << "\t" << color1 << endl;
//		cout << bub << endl;
		if(bubend[bub] == true){//meanning bub is a level-1 bubble
			bubfather[bub+1] = -1;
			continue;
		}
		for(int bub2 = bub;bub2 < bubnum;bub2++){
			if(bub2 == bub)
				continue;
			string color2 = bubcolors[bub2];
			if(color1 == color2)
				continue;
			vector<string> cmpcolorset = {color1,color2};
			string commoncolor = ColorCap(cmpcolorset);
			if(commoncolor == color1){
				bubfather[bub+1] = bub2+1;
//				cout << bub2+1 << "->" << bub+1 << endl;
				if(bubchild.find(bub2+1) != bubchild.end()){
					vector<int> tmpchild = bubchild[bub2+1];
					tmpchild.push_back(bub+1);
					bubchild[bub2+1] = tmpchild;
				}
				else{
					bubchild[bub2+1] = {bub+1};
				}
				break;
			}
//			if(color2 == stdcolor){
//				cout << "meet error during find fathers ~" << endl;
//			}
		}
	}
	cout << "Finished getting affiliation" << endl;
	map<int,int> bublevel; //store each bubble's nesting level
//	vector<int> bubvisited(bubnum,0);
	map<ssize_t,int> node2var;
	for(int bub = 0;bub < bubnum;bub++){//get each bubble's levels
		if(bubfather[bub+1] != -1)
			continue;
//		cout << "root bubble:" << bub << endl;

		int level = 1;
		vector<int> fatherbubid = {bub+1};
		while(!fatherbubid.empty()){
			int setnum = fatherbubid.size();
			vector<int> childidset;
			for(int bub1 = 0;bub1 < setnum;bub1++){
				int bubid = fatherbubid[bub1];
				bublevel[bubid] = level;
			//	cout << bubid << ":" << level << endl;
//				bubvisited[bubid-1] = 1;
				if(bubchild.find(bubid) != bubchild.end()){
					vector<int> childbubset = bubchild[bubid];
					childidset.insert(childidset.end(),childbubset.begin(),childbubset.end());
				}
			}
			fatherbubid = childidset;
			level++;
		}
	}
	cout << "Finished getting bubble level information" << endl;
	map<int,int> bub2vartype;//1:SNP;2:delete;3:Insert;4:Indel
	map<ssize_t,int> node2vartype;
	map<ssize_t,int> node2actualgap;
	int unsuregap = 0;
	int changecount = 0;
	string insertcolor = zerocolor;
	bool insertcond = false;
	int totalgaps = 0;
	for(int bub = 0;bub < bubnum;bub++){
		vector<ssize_t> pairnode = bubpairs[bub];
		ssize_t source = pairnode[0];
		if(node2vartype[source]){
			bub2vartype[bub+1] = node2vartype[source];
			continue;
		}
		string sourcecolor = ColorCap({getnodecolor(source),researched_samples});
		if(sourcecolor[numcolors - 1] != '1')
			continue;
		cout << "bub:" << bub+1 << "\t" << source << "\t" << insertcond <<  endl;
		if(!node2gap[source] ||node2gap[source] == 0){
			if(!node2vartype[source]){
				bub2vartype[bub+1] = 1;
				node2vartype[source] = 1;
			}else{
				bub2vartype[bub+1] = node2vartype[source];
			}
		}
		if(node2gap[source] > 0){
			string tmpgapcolor = node2gapcolor[source];
			int gaps = node2gap[source];
			if(insertcond == false){
				node2actualgap[source] = node2gap[source];
				if(tmpgapcolor[numcolors - 1] == '1'){
					node2vartype[source] = 3;
					bub2vartype[bub+1] = 3;
					insertcolor = sourcecolor;
					totalgaps = node2gap[source];
					insertcond = true;
					continue;
				}else{
					node2vartype[source] = 2;
					bub2vartype[bub+1] = 2;
				}
			}else{
				if(tmpgapcolor[numcolors - 1] == '1'){//insert
					node2vartype[source] = 3;
					bub2vartype[bub+1] = 3;
					totalgaps += gaps;
					//node2actualgap[source] = node2gap[source];//unsure
					node2actualgap[source] = totalgaps;//unsure
					insertcolor = ColorCup({insertcolor,sourcecolor});
				}else{//gapcolor doesn't include reference
					string capcolor = ColorCap({insertcolor,tmpgapcolor});
					cout << "totalgaps:" << totalgaps << "\t" << gaps << endl;
					if(capcolor == zerocolor){
						if(totalgaps < gaps){
							node2vartype[source] = 2;
							bub2vartype[bub+1] = 2;
							node2actualgap[source] = gaps - totalgaps;
						}else if(totalgaps == gaps){
							node2vartype[source] = 1;
							bub2vartype[bub+1] = 1;
							node2actualgap[source] = 0;
						}else{
							node2vartype[source] = 3;
							bub2vartype[bub+1] = 3;
							node2actualgap[source] = totalgaps - gaps;
						}
					}else{
						node2vartype[source] = 2;	
						bub2vartype[bub+1] = 2;
						node2actualgap[source] = gaps;
					}
					insertcolor = ColorCup({insertcolor,sourcecolor});
				}

			}
		
		}
		cout << "source:" << source << "\t" << node2vartype[source] << endl;
		if(insertcolor == researched_samples|| sourcecolor == researched_samples){
			insertcond = false;
			totalgaps = 0;
			insertcolor = zerocolor;
		}
	}
	cout << "Finished refsourcenodes' vartype detection" << endl;

	for(int bub = 0;bub < bubnum;bub++){
		vector<ssize_t> pairnode = bubpairs[bub];
		ssize_t source = pairnode[0];
		if(node2vartype[source] && node2vartype[source]!=0){
			bub2vartype[bub+1] = node2vartype[source];
			continue;
		}
//		cout << "bub:" << bub+1 << "\t" << source << "\t" << node2gap[source] << endl;
		if(!node2actualgap[source] ||node2actualgap[source] == 0){
			bub2vartype[bub+1] = 1;
			node2vartype[source] = 1;
		}
		if(node2gap[source] > 0){
//			cout << bub+1 << "\t" << source << "\t" << node2gap[source] << endl;
			bool getfather = false;
			int fathbub = bub+1;
			vector<ssize_t> gapnodes = {source};
			int fathsource = -1;
			while(getfather == false){//to find the nearest bubble of bub that includes reference or doesn't have gap
				int tmpbub = bubfather[fathbub];
				vector<ssize_t> pairnode2 = bubpairs[tmpbub-1];
				ssize_t source2 = pairnode2[0];
//					cout << "tmpbub:" << tmpbub << "\tsource2:" << source2 << endl;
				if(source2 != source){
					if(!node2vartype[source2]){
			//			string sourcecolor2 = ColorCap({getnodecolor(source2),researched_samples});
				//		if(node2gap[source2] > 0){//source2:indel
							if(std::find(gapnodes.begin(),gapnodes.end(),source2) == gapnodes.end()){//can't find the source2 in the gapnodes
								gapnodes.push_back(source2);
							}
				//		}else{//source2:snp
				//			getfather = true;
				//			fathsource = source2;
				//		}							
					}else{//father source's vartype exists
						getfather = true;
						fathsource = source2;
					}
				}
				fathbub = tmpbub;
			}
			int gapnodeslen = gapnodes.size();
//			if(gapnodeslen > 1)
//				cout << "fathersource:" << fathsource << "\tfathergapbub:" << fathbub << "\tdepth:" << gapnodeslen << "\trefgap:" << node2gap[fathsource] << endl;
			bool ref_ins = false;
			for(int i = gapnodeslen-1;i >= 0;i--){
				ssize_t childsource = gapnodes[i];
		//		int fathgap = node2gap[fathsource];
				int childgap = node2gap[childsource];
				int childvar = 1;
				int fathvar = 1;
		//		cout << fathsource << "->" << childsource << endl;
				if(node2vartype[fathsource]){//father source does't define the vartype
					fathvar = node2vartype[fathsource];
					if(childgap == 0){//1:SNP;2:Delete;3:Insert;4:Indel
						childvar = 1;
						node2actualgap[childsource] = 0;
						node2vartype[childsource] = childvar;
						fathsource = childsource;
						if(i == 0)
							bub2vartype[bub+1] = childvar;
					//	node2actualgap[childsource] = node2gap[childsource];
						continue;
						
					}
					if(i == gapnodeslen-1){
						if(fathvar <= 2){
							childvar = 2;
							node2actualgap[childsource] = node2gap[childsource];
						}else{
							childvar = 4;
							node2actualgap[childsource] = -node2gap[childsource];
							ref_ins = true;
							unsuregap++;
						}
					}else{
						if(fathvar <= 2){
							if(ref_ins == true){
								childvar = 4;
								node2actualgap[childsource] = -node2gap[childsource];
								unsuregap++;
							}else{
								childvar = 2;
								node2actualgap[childsource] = node2gap[childsource];
							}
						}else{
							childvar = 4;
							node2actualgap[childsource] = -node2gap[childsource];
							unsuregap++;
						}
					}
				}else{
					cout << "find father-child node vartype error!!!" << endl;
				}
				node2vartype[childsource] = childvar;
				fathsource = childsource;
				if(i == 0)
					bub2vartype[bub+1] = childvar;
			}
		}
	}
	cout << "Changed gap node count:" << changecount << endl;
	cout << "Unsure gap count:" << unsuregap << endl;;
	cout << "Finished getting simple vartype information" << endl;
//	map<int,vector<char>> pos2var;
//	map<int,int> pos2vartype;
///*	
	int range = maxnodepos + kmer-1;
	string zerostr(range,'-');
	vector<char> msamat(range*numcolors,'-');
	vector<string> msastr(numcolors,zerostr);
	for(ssize_t node = 0;node < numnodes;node++){
		if(!nodepos[node])
			continue;
		string nodelabel = dbg.node_label(node);
		string nodecolor = ColorCap({getnodecolor(node),researched_samples});
		if(nodepos[node] == 1){
			cout << "startnode:" << node << endl;
			nodecolor = researched_samples;
			for(int i = 0;i < numcolors;i++){
				if(nodecolor[i] == '1'){
					string sampstr = msastr[i];
					for(int j = 0;j < kmer;j++){
						sampstr[j] = nodelabel[j];
						msamat[i*range+j] = nodelabel[j];
					}
					msastr[i] = sampstr;
				}
			}
		}else{
			for(int i = 0;i < numcolors;i++){
				if(nodecolor[i] == '1'){
					string sampstr = msastr[i];
					sampstr[kmer - 2 + nodepos[node]] = nodelabel[kmer-1];
					msastr[i] = sampstr;
					msamat[i * range + kmer - 2 + nodepos[node]] = nodelabel[kmer-1];
				}
			}
		}
	}
	for(int i = 0;i < numcolors;i++){
		msafile << msastr[i] << endl;
	}
	cout << "Finished getting initial MSA result" << endl;
//*/
/*
	for(int i = kmer-1;i < range;i++){
		vector<char> posvar;
		for(int j = 0;j < numcolors;j++){
			cher tmpbase = msamat[j*range+i];
			if(std::find(posvar.begin(),posvar.end(),tmpbase) == posvar.end()){//can't find the tmpbase in the posvar
				posvar.push_back(tmpbase);
			}
		}
		pos2var[i+2-kmer] = posvar;
	}

*/
	map<int,vector<string>> refpos2var;
	map<int,int> refpos2nodepos;
	map<int,int> nodepos2refpos;
	map<int,int> refpos2vartype;
	int unsureref = 0;
	map<ssize_t,bool> node2visit;
	for(int bub = 0;bub < bubnum;bub++){
		vector<ssize_t> pairnode = bubpairs[bub];
		ssize_t source = pairnode[0];
		if(node2visit[source] == true)
			continue;
		ssize_t rootfathsource = -1;
		if(refnodepos[source] == -1){ //to find the nearest father source in reference
			bool getfather = false;
			int fathbub = bub+1;
			vector<ssize_t> gapnodes = {source};
			int fathsource = -1;
			while(getfather == false){//to find the nearest bubble of bub that includes reference 
				int tmpbub = bubfather[fathbub];
				vector<ssize_t> pairnode2 = bubpairs[tmpbub-1];
				ssize_t source2 = pairnode2[0];
//				cout << "tmpbub:" << tmpbub << "\tsource2:" << source2 << endl;
				if(source2 != source){
					//if(!source2used[source2]||refnodepos[source2] > -1)
					if(refnodepos[source2] != -1){
						getfather = true;
						fathsource = source2;
					}else{//does not find the reference bubble
						if(std::find(gapnodes.begin(),gapnodes.end(),source2) == gapnodes.end()){//can't find the source2 in the gapnodes
							gapnodes.push_back(source2);
						}
					}
				}
				fathbub = tmpbub;
			}
			int gapnodeslen = gapnodes.size();
			rootfathsource = fathsource;
//			if(gapnodeslen > 1)
//				cout << "fathersource:" << fathsource << "\tfathergapbub:" << fathbub << "\tdepth:" << gapnodeslen << "\trefgap:" << node2gap[fathsource] << endl;
			
//			cout << fathsource;
			for(int i = gapnodeslen-1;i >= 0;i--){
				ssize_t childsource = gapnodes[i];
//				cout << "->" << childsource;
				int fathgap = node2actualgap[fathsource];
				int childgap = node2actualgap[childsource];
				int fathvar = node2vartype[fathsource];
//				int childvar = node2vartype[childsource];
				int diff = nodepos[fathsource] - refnodepos[fathsource];
				if(fathvar <= 2){
					refnodepos[childsource] = nodepos[childsource] - diff;
				}else if(fathvar == 3){
					int gapdiff = fathgap - childgap;
					refnodepos[childsource] = nodepos[childsource] - diff - gapdiff;
				}else{
					if(fathgap > childgap){
						int gapdiff = fathgap - childgap;
						int rawpos = nodepos[childsource] - diff - gapdiff;
						refnodepos[childsource] = -rawpos;//- means not sure
						unsureref++;
					}else{
						int gapdiff = childgap - fathgap;
						int rawpos = nodepos[childsource] - diff + gapdiff;
						refnodepos[childsource] = -rawpos;
						unsureref++;
					}

				}
				nodepos2refpos[nodepos[childsource]] = abs(refnodepos[childsource]);
				refpos2nodepos[abs(refnodepos[childsource])] = nodepos[childsource];
				vector<string> posvar;
				int vartype = node2vartype[childsource];
				int gap = abs(node2actualgap[childsource]);
				int visit_range;
				if(!refpos2vartype[abs(refnodepos[childsource])+1]){
					refpos2vartype[abs(refnodepos[childsource])+1] = vartype;
				}else if(refpos2vartype[abs(refnodepos[childsource])+1] != vartype){
					//cout << "same position has different vartype1" << endl;
					cout << "same position has different vartype1:" << childsource << "\t" << nodepos[childsource] << "\t" << abs(refnodepos[childsource])+1 << "\t" << refpos2vartype[abs(refnodepos[childsource])+1] << "\t" << vartype << endl;
					refpos2vartype[abs(refnodepos[childsource])+1] = -(refpos2vartype[abs(refnodepos[childsource])+1] + vartype);
					//refpos2vartype[abs(refnodepos[childsource])+1] = refpos2vartype[abs(refnodepos[childsource])+1]> vartype? refpos2vartype[abs(refnodepos[childsource])+1]:vartype;
				}
				if(vartype == 1){
					visit_range = 1;
				}else if(vartype == 2){
					visit_range = gap;
				}else{
					visit_range = gap + 1;
				}
				vector<string> varstr(numcolors,"");
				vector<ssize_t> nodeset;
				queue<ssize_t> queueset;
				queueset.push(childsource);
				int startpos = nodepos[childsource];
		//		string sourcecolor = ColorCap({getnodecolor(childsource),researched_samples});
				//cout << "source:" << source << "\t" << startpos << endl;
				while(!queueset.empty()){
					ssize_t node = queueset.front();
					for(int k = 1;k < 5;k++){
						ssize_t nextnode = dbg.outgoing(node,k);
						if(nextnode == -1)
							continue;
						if(vartype > 1){
							if(nodepos[nextnode] - startpos > visit_range)
								continue;
						}
						string nextnodecolor = ColorCap({getnodecolor(nextnode),researched_samples});
						if(nextnodecolor == zerocolor)
							continue;
						string nodelabel = dbg.node_label(nextnode);
						for(int i = 0;i < numcolors;i++){
							if(nextnodecolor[i] == '1'){
								string tmpstr = varstr[i];
								tmpstr += nodelabel[kmer-1];
								varstr[i] = tmpstr;
						
							}
						}
						if(vartype > 1)
							queueset.push(nextnode);
					}
					queueset.pop();
				}
				string refstr;
				for(int i = 0;i < numcolors;i++){
					if(varstr[i] == "")
						continue;
					if(std::find(posvar.begin(),posvar.end(),varstr[i]) == posvar.end()){
						posvar.push_back(varstr[i]);
				//		cout << varstr[i] << endl;
					}
					if(i == numcolors - 1)
						refstr = varstr[i];
				}
//				refpos2var[refnodepos[childsource]] = posvar;
				int refsourcepos = abs(refnodepos[childsource])+1;
				if(!refpos2var[refsourcepos].empty()){
					vector<string> tmpposvar = refpos2var[refsourcepos];
				//	int strsize = tmpposvar.size();
			//		if(strsize == visit_range){
						int posvarsize = tmpposvar.size();
						for(int i = 0;i < posvarsize - 1;i++){
							if(std::find(posvar.begin(),posvar.end(),tmpposvar[i]) == posvar.end()){
								posvar.push_back(tmpposvar[i]);
							}
						}
			//		}
				}
				if(refstr.empty()){
					int tmprange = 1;
					if(vartype == 2){
						tmprange = gap;
					}
					refstr = refstring.substr(refsourcepos+kmer-1-1,tmprange);
					cout << refsourcepos+kmer-1-1 << "\t" << refstr<< endl;
				}
		//		cout << "refstr1:" << refstr << endl;
		//		if(std::find(posvar.begin(),posvar.end(),refstr) == posvar.end())
		//			posvar.push_back(refstr);
				refpos2var[refsourcepos] = posvar;
				refpos2var[refsourcepos].push_back(refstr);
				node2visit[childsource] = true;
				fathsource = childsource;
			}
			
		}
//		cout << endl;
		if(rootfathsource != -1 || refnodepos[source] != -1){
			if(rootfathsource != -1)
				source = rootfathsource;
	//		cout << source << "\t";
			if(node2visit[source] == true)
				continue;
			nodepos2refpos[nodepos[source]] = refnodepos[source];
			refpos2nodepos[refnodepos[source]] = nodepos[source];
			vector<string> posvar;
			int vartype = node2vartype[source];
	//		cout << vartype << endl;
			if(!refpos2vartype[abs(refnodepos[source])+1]){
				refpos2vartype[abs(refnodepos[source])+1] = vartype;
			}else if(refpos2vartype[abs(refnodepos[source])+1] != vartype){
				cout << "same position has different vartype2:" << source << "\t" << nodepos[source] << "\t" << abs(refnodepos[source])+1 << "\t" << refpos2vartype[abs(refnodepos[source])+1] << "\t" << vartype << endl;
				refpos2vartype[abs(refnodepos[source])+1] = -(refpos2vartype[abs(refnodepos[source])+1] + vartype);
			//	refpos2vartype[abs(refnodepos[source])+1] = refpos2vartype[abs(refnodepos[source])+1]> vartype? refpos2vartype[abs(refnodepos[source])+1]:vartype;
			}
			int gap = abs(node2actualgap[source]);
		//	cout << "gap:" << gap << endl;
			int visit_range;
			if(vartype == 1){
				visit_range = 1;
			}else if(vartype == 2){
				visit_range = gap;
			}else{
				visit_range = gap + 1;
			}
		//	if(source == 38146){
		//		cout << source << "\t" << vartype << "\t" << gap << "\t" << visit_range << endl;
		//		cout << node2gapcolor[source] << endl;
		//	}
			vector<string> varstr(numcolors,"");
			vector<ssize_t> nodeset;
			queue<ssize_t> queueset;
			queueset.push(source);
			int startpos = nodepos[source];
			while(!queueset.empty()){
				ssize_t node = queueset.front();
				for(int k = 1;k < 5;k++){
					ssize_t nextnode = dbg.outgoing(node,k);
					if(nextnode == -1)
						continue;
					if(vartype > 1){
						if(nodepos[nextnode] - startpos > visit_range)
							continue;
					}
					string nextnodecolor = ColorCap({getnodecolor(nextnode),researched_samples});
					if(nextnodecolor == zerocolor)
						continue;
					string nodelabel = dbg.node_label(nextnode);
					for(int i = 0;i < numcolors;i++){
						if(nextnodecolor[i] == '1'){
							string tmpstr = varstr[i];
							tmpstr += nodelabel[kmer-1];
							varstr[i] = tmpstr;
					
						}
					}
					if(vartype > 1)
						queueset.push(nextnode);
				}
				queueset.pop();
			}
			string refstr;
			for(int i = 0;i < numcolors;i++){
				if(varstr[i] == "")
					continue;
				if(std::find(posvar.begin(),posvar.end(),varstr[i]) == posvar.end()){
					posvar.push_back(varstr[i]);
				}
				if(i == numcolors - 1)
					refstr = varstr[i];
			}
		//	if(source == 38146){
		//		cout << refstr << endl;
		//	}
			int refsourcepos = abs(refnodepos[source])+1;
			if(!refpos2var[refsourcepos].empty()){
				vector<string> tmpposvar = refpos2var[refsourcepos];
				int posvarsize = tmpposvar.size();
				for(int i = 0;i < posvarsize - 1;i++){
					if(std::find(posvar.begin(),posvar.end(),tmpposvar[i]) == posvar.end()){
						posvar.push_back(tmpposvar[i]);
					}
				}
			}
			if(refstr.empty()){
				int tmprange = 1;
				if(vartype == 2){
					tmprange = gap;
				}
				refstr = refstring.substr(refsourcepos+kmer-1-1,tmprange);
					cout << refsourcepos+kmer-1-1 << "->" << refstr<< endl;
			}
		//	cout << "refstr2:" << refstr << endl;
		//	if(std::find(posvar.begin(),posvar.end(),refstr) == posvar.end())
		//		posvar.push_back(refstr);
			refpos2var[refsourcepos] = posvar;
			refpos2var[refsourcepos].push_back(refstr);
		}
		node2visit[source] = true;
	}
	cout << "Unsure refpos count:" << unsureref << endl;
	cout << "Finished getting variation result" << endl;

	varfile << "refpos\tnodepos\tref\talt\ttype" << endl;
	for(int pos = 1;pos < maxnodepos;pos++){//refpos
		if(!refpos2var[pos].empty()){
		//	cout << pos << endl;
			vector<string> vars = refpos2var[pos];
			int len = vars.size();
			int nodepos1 = refpos2nodepos[pos-1]+1;
			varfile << pos+kmer-1 << "\t" << nodepos1 << "\t" << vars[len-1] << "\t" << vars[0];
			for(int i = 1;i < len-1;i++){
				varfile << ":" << vars[i];
			}
			varfile << "\t" << refpos2vartype[pos];
			varfile << endl;
		}
	}

	map<int,int> bub2var;
	map<int,int> bub2gap;
	for(int bub = 0;bub < bubnum;bub++){ ///output and store the level/pairs/color/affiliation information
		vector<ssize_t> pairnodeset;
		vector<int> affiliationset;
		ssize_t level = bublevel[bub+1];
		pairnodeset.push_back(level);
		pairnodeset.insert(pairnodeset.end(),bubpairs[bub].begin(),bubpairs[bub].end());
		pairnodesets.push_back(pairnodeset);
		affiliationset.push_back(bub+1);
		affiliationset.push_back(bubfather[bub+1]);
		bubblecolors.push_back(bubcolors[bub]);
//		cout << bub << "\t" << bublevel[bub] << "\t";
		vector<ssize_t> tmppair = bubpairs[bub];
//		cout << tmppair[0] << ":" << tmppair[1] << "\t" << bubcolors[bub] << "\t" << bubfather[bub] << "\t";
		pairsfile << bub+1 << "\t" << level << "\t" << nodepos[tmppair[0]] << "\t" << nodepos[tmppair[1]] << "\t" << tmppair[0] << ":" << tmppair[1] << "\t";
		pairsfile << refnodepos[tmppair[0]] << "\t" << refnodepos[tmppair[1]] << "\t" << node2actualgap[tmppair[0]] << "\t" << bub2vartype[bub+1]<< "\t";
		pairsfile << bubfather[bub+1] << ";";
		bub2var[bub+1] = bub2vartype[bub+1];
		bub2gap[bub+1] = node2actualgap[tmppair[0]];
//		affiliationfile << bub << ";" << bubfather[bub] << ";"; 
		if(bubchild.find(bub+1) != bubchild.end()){
			vector<int> childidset = bubchild[bub+1];
			affiliationset.insert(affiliationset.end(),childidset.begin(),childidset.end());
			int childnum = childidset.size();
//			affiliationfile << childidset[0];
			pairsfile << childidset[0];
//			cout << childidset[0];
			for(int i = 1;i < childnum;i++){
//				cout << ":" << childidset[i];
//				affiliationfile << ";" << childidset[i];
				pairsfile << ";" << childidset[i];
			}
//			affiliationfile << endl;
			pairsfile << "\t";
		}
		else{
			affiliationset.push_back(-1);
//			cout << -1;
//			affiliationfile << -1 << endl;
			pairsfile << -1 << "\t";
		}
		bool getbootfather = false;
		string bubstructure = to_string(bub+1);
		int tmpbub = bub+1;
		while(getbootfather == false){
			int bubfath = bubfather[tmpbub];
			if(bubfath > 0){
				string tmpstr = to_string(bubfath);
				bubstructure = "_" + bubstructure;
				bubstructure = tmpstr + bubstructure;
				tmpbub = bubfath;
			}else{
				bubstructure = "bub" + bubstructure;
				getbootfather = true;
			}
		}
		pairsfile << bubstructure << "\t";
		pairsfile << bubcolors[bub] << endl;;
//		cout << endl;
		affiliationsets.push_back(affiliationset);
	}
//	int num_nodes = dbg.num_nodes();
//	for(int startnode = 0;startnode < num_nodes;startnode++){
//		cout << startnode << "\t" << nestlevelsets[startnode] << endl;
//		nestlevelfile << startnode << "\t" << nestlevelsets[startnode] << endl;
//	}
	cout << "Finished output mypairnodes work" << endl;

	colored_su.pairnodesets = pairnodesets;
	colored_su.affiliationsets = affiliationsets;
	colored_su.bubblecolors = bubblecolors;
	colored_su.nodepos = nodepos;
	colored_su.refnodepos = refnodepos;
	colored_su.bub2gap = bub2gap;
	colored_su.bub2var = bub2var;
//	colored_su.researched_samples = researched_samples;
	int totaltime = getMilliSpan(t);
	if(totaltime < 1000){
		cout << "Total calculation time: " << getMilliSpan(t) << "ms" << endl;
	}
	else{
		int tolsecond = totaltime/1000;
		int hour = tolsecond/3600;
		int second1 = tolsecond%3600;
		int minute = second1/60;
		int second = second1%60;
		cout << "Total get colored-superbubble calculation time: " << hour << "h" << minute << "m" << second << "s" << endl;
	}
	return colored_su;
}

//Get the nestinglevelpairs from the nestinglevelpairsfile  
//File format:level	node1:node2	color fatherbub:child1:child2:child3
C_Su read_mypairs_from_file(string csupb_topo_file){
	C_Su c_su;
	ifstream pairsfile; 
    pairsfile.open(csupb_topo_file.data());    
	ssize_t level;
	string pairnodes,colors,affinfo,bubstructure;
	vector<vector<ssize_t>> csupb_topo_set;
	vector<vector<int>> affiliationsets;
	vector<string> bubblecolors;
    if(!pairsfile.is_open()){ //if read fail,break
		cout << "open csupb_topo_file failure!" << endl;
	}
	else{
		cout << "open pairs file successfully!" << endl;
	}
	int order = 1;
	int bubid;
	int startpos;
	int endpos;
	int refstartpos,refendpos,gaps,var;
	map<ssize_t,int> nodepos;
	map<ssize_t,int> refnodepos;
	map<int,int> bub2gap;
	map<int,int> bub2var;
//	map<ssize_t,int> nodepos = read_mynodepos_from_file("mynodeposinfo.txt");
	while(pairsfile >> bubid >> level >> startpos >> endpos >> pairnodes >> refstartpos >> refendpos >> gaps >> var >> affinfo >> bubstructure >> colors){
//		cout << bubid << "\t" << level << "\t" << startpos << "\t" << endpos << "\t" << pairnodes << "\t" << affinfo << "\t" << bubstructure << "\t" << colors << endl;
		int strlen = pairnodes.size();
		bub2gap[bubid] = gaps;
		bub2var[bubid] = var;
		vector<int> nodes;
		int start = 0;
		for(int i = 0;i < strlen;i++){
			if(pairnodes[i] == ':'){
				int idlen = i-start;
				string id = pairnodes.substr(start,idlen);
				int intid = atoi(id.c_str());
				nodes.push_back(intid);
				//ssize_t sourcenode = intid;
	//			nodepos[sourcenode] = startpos;
				start = i + 1;
			}
			if(i == strlen - 1){
				int idlen = i-start+1;
				string id = pairnodes.substr(start,idlen);
				int intid = atoi(id.c_str());
			//	ssize_t sinknode = intid;
	//			nodepos[sinknode] = endpos;
				nodes.push_back(intid);
			}
		}
		int lennode = nodes.size();
//		cout<< "nodes counts:" << lennode << endl;
		vector<ssize_t> csupbpair;
		csupbpair.push_back(level);
		for(int m = 0;m < lennode;m++){
			ssize_t nodesid = nodes[m];
			csupbpair.push_back(nodesid);
		}
		csupb_topo_set.push_back(csupbpair);
		bubblecolors.push_back(colors);
		int affstrlen = affinfo.size();
		vector<int> info;
		info.push_back(order);
		order++;
		int affstart = 0;
		for(int i = 0;i < affstrlen;i++){
			if(affinfo[i] == ';'){
				int idlen = i-affstart;
				string id = affinfo.substr(affstart,idlen);
				int intid = atoi(id.c_str());
				info.push_back(intid);
				affstart = i + 1;
			}
			if(i == affstrlen - 1){
				int idlen = i-affstart+1;
				string id = affinfo.substr(affstart,idlen);
				int intid = atoi(id.c_str());
				info.push_back(intid);
			}
		}
		affiliationsets.push_back(info);
	}
    pairsfile.close();
	c_su.pairnodesets = csupb_topo_set;
	c_su.affiliationsets = affiliationsets;
	c_su.bubblecolors = bubblecolors;
	c_su.bub2gap = bub2gap;
	c_su.bub2var = bub2var;
	ifstream nodeposfile;
	string nodeposfilename = "offsetinfo.txt";
	nodeposfile.open(nodeposfilename.data());
	ssize_t mynode;
	int pos,refpos,gap;
	if(!nodeposfile.is_open()){
		cout << "open offsetinfo failure!" << endl;
	}
	else{
		while(nodeposfile >> mynode >> pos >> refpos >> gap){
			if(nodepos[mynode] && nodepos[mynode] < pos){
				cout << "offset value file has error" << endl;
			}
			nodepos[mynode] = pos;
			refnodepos[mynode] = refpos;
		//	if(mynode == 33745 || mynode == 9530){cout << mynode << ":" << pos << endl;}
		}
	}
	c_su.nodepos = nodepos;
	c_su.refnodepos = refnodepos;
	cout << "Finish loading data from file(mypairfile/mynodeposfile)" << endl;
	return c_su;
}


//#####################################################################################################################################
//BubInfo BFSget_Bub_Info(vector<int> nestlevelset,BubInfo bubbleinfo){//,vector<ssize_t> pairnodelist){
BubInfo BFSget_Bub_Info(BubInfo bubbleinfo){//,vector<ssize_t> pairnodelist){
//	cout << "bubid:" << bubbleinfo.bubbleid << endl;
//	ssize_t node = bubbleinfo.firstnode;
	vector<ssize_t> pairnodeset = bubbleinfo.pairnodeset;
//	cout << pairnodeset[1] << ":" << pairnodeset[2] << endl;
	map<ssize_t,int> nodepos = bubbleinfo.nodespos;
	ssize_t sourcenode = pairnodeset[1];
	ssize_t sinknode = pairnodeset[2];
	int sourcepos = nodepos[sourcenode];
	int sinkpos = nodepos[sinknode];
//	cout << pairnodeset[1] << ":" << pairnodeset[2] << endl;
//	string varpos(sinkpos - sourcepos,'0');
	string bubcolor = bubbleinfo.bubcolor;
//	string researched_samples = bubbleinfo.researched_samples;
//	cout << sourcenode << ":" << sinknode << "\t" << sourcepos << ":" << sinkpos << endl;
	bool isstart = true;
	queue<ssize_t> nodelist;
	nodelist.push(sourcenode);
	string nodelabel = dbg.node_label(sourcenode);
	int kmer = dbg.k-1;
	vector<int> visited(dbg.num_nodes(),0);
	visited[sourcenode] = 1;//0:unvisited 1:visiting 2:visited 3:incoming node unvisited completely
	//string tmpcolor1 = getnodecolor(sourcenode);
	//int numcolors = tmpcolor1.size();
	int numcolors = colors.size()/dbg.size();
	vector<int> emptyvec;
	string emptystr = "";
	vector<int> sampleweight(numcolors,0);
	vector<int> samplelength(numcolors,0);
//	vector<int> samplevirtuallength(numcolors,0);
	map<int,int> basisgap;//basisid->gaplen
	vector<vector<int>> samplebasis(numcolors,emptyvec);
	vector<Basis> bubblebasis;
//	vector<string> variantpos(numcolors,emptystr);
	string stdcolor(numcolors,'1');
	string zerocolor(numcolors,'0');

	vector<string> pathseqs(numcolors,"1");
	vector<string> pathcolors(numcolors,"1");
	string tmpstartstr = nodelabel.substr(0,kmer-1);
	for(int i = 0;i < numcolors;i++){
		if(bubcolor[i] == '1'){
			pathseqs[i] = tmpstartstr;
			pathcolors[i] = bubcolor;
		}
	}
	vector<int> pathlengths(numcolors,0);
	vector<int> pathsupernodes(numcolors,0);
	vector<int> pathweights(numcolors,0);
	vector<vector<int>> pathbranchwidths(numcolors,emptyvec);
//	vector<ssize_t> pathendnodes;
//	map<ssize_t,int> pathsourcenodes;
	vector<ssize_t> sourcenodes;
	vector<ssize_t> sinknodes;
	vector<vector<int>> branchdigits;
//	vector<int> samplediff(numcolors*numcolors,0);
	vector<vector<int>> weightmat(numcolors,emptyvec);
//	cout << "start to visited colored-superbubble " << node << endl;
	int basisid = 0;
//	vector<int> gapbasisid;
	while(!nodelist.empty()){
		ssize_t startnode = nodelist.front();
		int indegree = dbg.indegree(startnode);
		if(indegree > 1)
			sinknodes.push_back(startnode);
		int outdegree = dbg.outdegree(startnode);
		if(outdegree > 1)
			sourcenodes.push_back(startnode);
//		cout << "startnode:" << startnode << endl;
	//	visited[startnode]s = 1;
		int value;
		string color = getnodecolor(startnode);
		string commoncolor = ColorCap({color,bubcolor});
		Basis basisinfo;
		basisinfo.type = 1;
		basisinfo.startnode = startnode;
		basisinfo.endnode = startnode;
		basisinfo.color = commoncolor;
		basisinfo.length = 1;
		vector<ssize_t> incomingnodes;
		vector<ssize_t> outgoingnodes;
		for(unsigned long x = 1;x < dbg.sigma + 6;x++){
			ssize_t neighbornode;
			if(x < 5){
				neighbornode = dbg.outgoing(startnode,x);
				if(neighbornode == -1)
					continue;
				outgoingnodes.push_back(neighbornode);
			}
			else{
				neighbornode = dbg.incoming(startnode,x-5);
				if(neighbornode == -1)
					continue;
				incomingnodes.push_back(neighbornode);
			}
		}
		basisinfo.incomingnode = incomingnodes;
		basisinfo.outgoingnode = outgoingnodes;
		string supernodelabel = dbg.node_label(startnode);
		basisinfo.sequence = supernodelabel;
		int weight;
		int basiswidth = 1;
//		vector<int> virtuallen;
//		if(commoncolor == stdcolor){
//			//weight = kmer * numcolors;
//			weight = kmer * numcolors;
//			basiswidth = numcolors;
//		}
//		else{
			vector<int> coloredsamp;
			vector<int> uncoloredsamp;
			for(int i = 0;i < numcolors;i++){
				if(commoncolor[i] == '1'){
					coloredsamp.push_back(i);
		//			if(samplevirtuallength[i] > 0)
		//				virtuallen.push_back(samplevirtuallength[i]);
				}
				else{uncoloredsamp.push_back(i);}
			}
			weight = kmer * coloredsamp.size();	
			basiswidth = coloredsamp.size();
			//int uncolorednum = uncoloredsamp.size();	
//		}
//		int maxvalue = Max(virtuallen);
//		int minvalue = Min(virtuallen);
		basisinfo.width = basiswidth;
		bubblebasis.push_back(basisinfo);
//		bool havegap = false;
//		cout << "virtuallen:" << virtuallen << endl;
//		cout << "maxvalue:" << maxvalue << "\tminvalue:" << minvalue << endl;
//		if(maxvalue > minvalue && minvalue > 0){
//			havegap = true;
//			for(int i = 0;i < basiswidth;i++){
//				int sampid = coloredsamp[i];
//				if(samplevirtuallength[sampid] < maxvalue){
//					int diffvalue = maxvalue-samplevirtuallength[sampid];
//					samplevirtuallength[sampid]+= diffvalue;
//					vector<int> sampbasis = samplebasis[sampid];
//					int lastbasisid = sampbasis[sampbasis.size()-1];
//					basisgap[lastbasisid] = diffvalue;
				//			auto it2 = std::find(std::begin(boundnodeset)+1,std::end(boundnodeset),pos);
				//			if(it2 == std::end(boundnodeset)){//can not find node
//					auto it1 = std::find(std::begin(gapbasisid),std::end(gapbasisid),lastbasisid);
//					if(it1 == std::end(gapbasisid))
//						gapbasisid.push_back(lastbasisid);
//					string varpos = variantpos[sampid];
//					for(int j = 0;j < diffvalue;j++){
//						varpos += '2';
//					}
//					variantpos[sampid] = varpos;
//				}
//			}
//		}
		for(int i = 0;i < numcolors;i++){
			vector<int> weightvec = weightmat[i];
			if(commoncolor[i] == '1'){
				weightvec.push_back(weight);
				samplelength[i]++;
				//samplevirtuallength[i]++;
				sampleweight[i]+=basiswidth;
				vector<int> sampbasis = samplebasis[i];
				sampbasis.push_back(basisid);
				samplebasis[i] = sampbasis;
	//			string varpos = variantpos[i];
	//			if(dbg.outdegree(startnode) > 1 && startnode != sinknode)
				//if(havegap == true)
	//				varpos += '1';
	//			else{varpos += '0';}
	//			variantpos[i] = varpos;
			
				pathweights[i] += weight;
				pathlengths[i]++;
				pathcolors[i] = ColorCap({pathcolors[i],commoncolor});
				pathsupernodes[i]++;
				pathseqs[i] += supernodelabel.substr(kmer-1,1);

				pathbranchwidths[i].push_back(basiswidth);
			}
			else{weightvec.push_back(0);}
			weightmat[i] = weightvec;
		}
		basisid++;

		vector<int> outnodepos;
//		int startnodepos = nodepos[startnode];
		for(unsigned long x = 1;x < dbg.sigma + 1;x++){
			if(startnode == sinknode){
			//	branchdigits.push_back(branchdigit);
				break;
			}
			vector<int> branchdigit;
			string boundnodestr = dbg.node_label(startnode);
			int strlen = boundnodestr.size();
			for(int k = 0;k < strlen;k++){
				char element = boundnodestr[k];
				if(element == 'A')
					branchdigit.push_back(1);
				else if(element == 'C')
					branchdigit.push_back(2);
				else if(element == 'G')
					branchdigit.push_back(3);
				else if(element == 'T')
					branchdigit.push_back(4);
				else{branchdigit.push_back(0);}
			}
			ssize_t pos;
			pos = dbg.outgoing(startnode,x);
			if(pos == -1)
				continue;
			if(visited[pos] == 2)
				continue;
			string poscolor = getnodecolor(pos);
			string capcolor = ColorCap({poscolor,bubcolor});
			if(capcolor == zerocolor)
				continue;
			int tmpnodepos = nodepos[pos];
			outnodepos.push_back(tmpnodepos);
			bubbleinfo.numedges++;
//			cout << node << "=>" << pos <<endl;
			string branch_labels = dbg.node_label(startnode);
			vector<ssize_t> branchnode;
		//	string color = getnodecolor(startnode);
			vector<string> branchcolor;
			value = x;
			branchnode.push_back(startnode);
			branchcolor.push_back(commoncolor);
			branch_labels += base[x];
			branchnode.push_back(pos);
			branchcolor.push_back(capcolor);
			vector<ssize_t> tmplist;
			if(visited[pos] == 0 && isstart == false){
				tmplist = {pos};
			}
			else if(visited[pos] == 0 && isstart == true){
				tmplist = {startnode,pos};
				isstart = false;
			}
			int nodenum = tmplist.size();
			for(int k = 0;k < nodenum;k++){
				pos = tmplist[k];
				color = ColorCap({bubcolor,getnodecolor(pos)});
				//cout << "pos: " << pos << " color:" << color << endl;
				if(k == 0)
					bubbleinfo.numnodes++;
			}
			int branchlen = 0;
			ssize_t branchstart = pos;
			ssize_t branchend = pos;
			string branchcol;
			bool normalnode = false;
			if(myindegree(pos,bubcolor) == 1 && myoutdegree(pos,bubcolor) == 1)
				normalnode = true;
			while(normalnode == true){
//			cout << "pos" << pos << "\tits degree:" << dbg.indegree(pos) << "\t" << dbg.outdegree(pos) << " " << bubbleinfo.visited[pos] << endl;
				branchlen++;
				branchend = pos;
				if(branchlen == 1)
					branchcol = ColorCap({getnodecolor(pos),bubcolor});
				visited[pos] = 2;
				ssize_t nextnode;
				branchdigit.push_back(value);
				for(unsigned long x2 = 1;x2 < dbg.sigma + 1;x2++){
			    	nextnode = dbg.outgoing(pos, x2);
					if(nextnode == -1)
						continue;
					string nextnodecolor = getnodecolor(nextnode);
					string tmpcapcolor = ColorCap({bubcolor,nextnodecolor});
					if(tmpcapcolor == zerocolor)
						continue;
					bubbleinfo.numedges++;
					pos = nextnode;
					value = x2;
//				cout << x2 << "\t" << pos << "\t" << nextnode << "\t" << visited[nextnode] << endl;
					branch_labels = branch_labels + base[x2];
					branchnode.push_back(pos);
				//	color = getnodecolor(pos);
					branchcolor.push_back(tmpcapcolor);
					if(visited[nextnode] == 0 && visited[nextnode] != 3){
						bubbleinfo.numnodes++;
						color = tmpcapcolor;
					}
					if(!((dbg.indegree(pos) == 1 && dbg.outdegree(pos)==1)||(myindegree(pos,bubcolor)==1 && myoutdegree(pos,bubcolor)==1)))
						normalnode = false;
					break;
				}
		//		bubbleinfo.visited[pos] = 1;
			}
			branchdigit.push_back(value);
			if(branchlen > 0){
				Basis branchbasisinfo;
				branchbasisinfo.type = 2;
				branchbasisinfo.startnode = branchstart;
				branchbasisinfo.endnode = branchend;
				vector<ssize_t> branchincomingnode = {startnode};
				vector<ssize_t> branchoutgoingnode = {pos};
				branchbasisinfo.incomingnode = branchincomingnode;
				branchbasisinfo.outgoingnode = branchoutgoingnode;
				branchbasisinfo.color = branchcol;
				branchbasisinfo.length = branchlen;
				branchbasisinfo.sequence = branch_labels;
				vector<int> coloredsamp1;
				for(int i = 0;i < numcolors;i++){
					if(branchcol[i] == '1')
						coloredsamp1.push_back(i);
				}
//				weight = (kmer - 1 + branchlen) * coloredsamp1.size();	
				weight = branchlen*kmer * coloredsamp1.size();	
//				weight = (branchlen+kmer-1) * coloredsamp1.size();	
				basiswidth = coloredsamp1.size();
				branchbasisinfo.width = basiswidth;
				bubblebasis.push_back(branchbasisinfo);
				for(int i = 0;i < numcolors;i++){
					vector<int> weightvec = weightmat[i];
					if(branchcol[i] == '1'){
						weightvec.push_back(weight);
						//sampleweight[i] += basiswidth;
						sampleweight[i] += weight;
						samplelength[i]+= branchlen;
						vector<int> sampbasis = samplebasis[i];
						sampbasis.push_back(basisid);
						samplebasis[i] = sampbasis;

						pathcolors[i] = ColorCap({pathcolors[i],branchcol});
						pathlengths[i] += branchlen;
						pathweights[i] += weight;
				//		pathendnodes[i] = pos;
						pathseqs[i] += branch_labels.substr(kmer,branch_labels.size() - kmer - 1);
						pathbranchwidths[i].push_back(basiswidth);
					}
					else{weightvec.push_back(0);}
					weightmat[i] = weightvec;
				}
				basisid++;
			}
		
//		cout << "supernode:" << pos << endl;
			bubbleinfo.numbranches++;
			bubbleinfo.branches.push_back(branch_labels);
			bubbleinfo.branchnodes.push_back(branchnode);
			bubbleinfo.branchcolors.push_back(branchcolor);
			branchdigits.push_back(branchdigit);
			branchcol = branchcolor[1];
			int branchsize = branchcolor.size();
			if(branchsize <= 2){
				for(unsigned long x = 0;x < 5;x++){
					ssize_t prenode = dbg.incoming(branchnode[1],x);
					if(prenode == -1)
						continue;
					if(prenode != branchnode[0])
						continue;
					ssize_t edge = dbg.incoming_edge(branchnode[1],x);
					branchcol = getedgecolor(edge);
					branchcol = ColorCap({branchcol,researched_samples});
					break;
				}
			}
			bool innodevisited = true;
			for(unsigned long x2 = 0;x2 < dbg.sigma+1;x2++){
				ssize_t frontnode = dbg.incoming(pos,x2);
				if(frontnode == -1)
					continue;
				string branchcol = ColorCap({getnodecolor(frontnode),bubcolor});
				if(branchcol == zerocolor || frontnode == startnode)
					continue;
				if(visited[frontnode] != 2){
					innodevisited = false;
					visited[pos] = 3;

					break;
				}
			}
			if(innodevisited == true){//means all the incoming nodes of sink node have been visited
//			if(visited[pos] == 0){
				bubbleinfo.numsupernodes++;
//			cout << "supernode:" << pos << endl;
//				bubbleinfo.firstnode = pos;
				nodelist.push(pos);
				visited[pos] = 1;
			}
		}
		
		visited[startnode] = 2;
		nodelist.pop();
/*
		int outnodepossize = outnodepos.size();
		int minpos = nodepos[sinknode];
		int maxpos = nodepos[sourcenode];
		if(outnodepossize > 0){
			for(int i = 0;i < outnodepossize;i++){
				int tmpnodepos = outnodepos[i];
				if(tmpnodepos < minpos)
					minpos = tmpnodepos;
				if(tmpnodepos > maxpos)
					maxpos = tmpnodepos;
			}
			if(minpos > startnodepos + 1){
				int startvarpos = startnodepos - sourcepos;
				for(int i = startvarpos;i < minpos - sourcepos - 1;i++){
					if(varpos[i] == '0')
						varpos[i] = '2';
					else if(varpos[i] == '1')
						varpos[i] = '3';
				}
			}
			if(minpos < startnodepos + 1){
				cout << "error when getting  the max/min position" << endl;
				cerr << "error when getting  the max/min position" << endl;
			}
			if(minpos == maxpos){
				if(outnodepossize > 1){
					int tmpvarpos = minpos - sourcepos - 1;
					if(varpos[tmpvarpos] == '0')
						varpos[tmpvarpos] = '1';
					else if(varpos[tmpvarpos] == '2')
						varpos[tmpvarpos] = '3';
				}
			}else{
				for(int i = 0;i < maxpos-minpos;i++){
					int tmpvarpos = minpos - sourcepos - 1 + i;
					if(varpos[tmpvarpos] == '0')
						varpos[tmpvarpos] = '2';
					else if(varpos[tmpvarpos] == '1')
						varpos[tmpvarpos] = '3';
				}
				if(maxpos - minpos > 50){
					cerr << "meet SV at the node " << sourcepos << endl;
					cout << "meet SV \nsource:sink" << sourcenode << ":" << sinknode << endl;
				}
			}
		}
*/				
	}
/*
	string::iterator it;
	it = varpos.end()-1;
	bool delgetend = false;
	while(delgetend == false){
		if(*it == '0'){
			varpos.erase(it);
			it--;
		}else{
			delgetend = true;
		}
	}
*/
	bubbleinfo.bubblebasis = bubblebasis;
	bubbleinfo.weightmat = weightmat;
	bubbleinfo.sourcenodes = sourcenodes;
	bubbleinfo.sinknodes = sinknodes;
	bubbleinfo.branchdigits = branchdigits;
//	bubbleinfo.variantspos = varpos;
	bubbleinfo.sampleweight = sampleweight;
	bubbleinfo.samplelength = samplelength;
	bubbleinfo.maxlength = sinkpos - sourcepos + 1;
//	bubbleinfo.basiswidths = basiswidths;
	bool unique = false;
	int uniquepos = 0;
	cout << "Start to delete repeat elements" << endl;
	while(unique == false){
		int pathnum = pathseqs.size();
		if(uniquepos >= pathnum - 1){
			unique = true;
			continue;
		}
		for(int i = pathnum - 1;i > uniquepos;i--){
			if(pathseqs[i] == pathseqs[uniquepos] || pathweights[i] == 0){
				pathseqs.erase(pathseqs.begin() + i);
				pathweights.erase(pathweights.begin() + i);
				pathcolors.erase(pathcolors.begin() + i);
				pathsupernodes.erase(pathsupernodes.begin() + i);
				pathlengths.erase(pathlengths.begin() + i);
				pathbranchwidths.erase(pathbranchwidths.begin() + i);
			}
		}
		uniquepos++;
	}
	if(pathweights[0] == 0){
		pathseqs.erase(pathseqs.begin());
		pathweights.erase(pathweights.begin());
		pathcolors.erase(pathcolors.begin());
		pathsupernodes.erase(pathsupernodes.begin());
		pathlengths.erase(pathlengths.begin());
		pathbranchwidths.erase(pathbranchwidths.begin());
	}
//	int pathnum = pathseqs.size();
//	cout << "unique sequence number:" << pathnum << endl;
//	for(int i = 0;i < pathnum;i++){
//		cout << pathseqs[i] << "\t" << pathcolors[i] << "\t" << pathlengths[i] << "\t" << pathsupernodes[i] << "\t" << pathweights[i] << endl;
//	}
	bubbleinfo.pathseqs = pathseqs;
	bubbleinfo.pathcolors = pathcolors;
	bubbleinfo.pathlengths = pathlengths;
	bubbleinfo.pathweights = pathweights;
	bubbleinfo.pathsupernodes = pathsupernodes;
	bubbleinfo.pathbranchwidths = pathbranchwidths;	
	return bubbleinfo;
}

vector<BubInfo> SupperBubblesAnalysis(C_Su c_su){//vector<int> nestlevelset,vector<vector<ssize_t>> pairnodeset){
	int t = getMilliCount();
//	vector<int> nestlevelset = c_su.nestlevelsets;
	vector<vector<ssize_t>> pairnodeset = c_su.pairnodesets;
	vector<string> bubblecolors = c_su.bubblecolors;
	vector<vector<int>> affinfo = c_su.affiliationsets;
	map<ssize_t,int> nodepos = c_su.nodepos;
//	string researched_samples = c_su.researched_samples;
	ofstream SupBubInfo("superbubbleinfo.txt");
//	ofstream AllsampInfo("allsampsuperbubble.txt");
	ofstream weightinfo("weightinfo.txt");
	ofstream basisinfo("basisinfo.txt");
	ofstream Branchwidth("pathbranchwidth.txt");
	ofstream simplesnpinfo("simplesnpinfo.txt");
	ofstream sampleweightinfo("sampleweight.txt");
	basisinfo << "bub\tbasisid\ttype\twidth\tlength\tstartnode\tendnode\tincomingnodes\toutgoingnodes\tsequence" << endl;
	ofstream Pathinfo("pathinfo.txt");
	Pathinfo << "bubbleid\tpathid\tpathwidth\tpathlen\tpathsupernodes\tactualweight\tpathweight\tpathseqs\tpathcolors" << endl;
	//bubbleid;pathid;pathwidth;pathlength;pathsupernodes;pathweights;path sequence;pathcolors
// samplelist can be used to classificate the node for special sample.
	vector<BubInfo> bubbleinfoset;
//	ssize_t num_nodes = dbg.num_nodes();
	int numcolors = colors.size()/dbg.size();
	string zerocolor(numcolors,'0');
	string stdcolor(numcolors,'1');
	int setnum = pairnodeset.size();
	vector<int> emptyvec;
	vector<vector<int>> snpvec(numcolors,emptyvec);
/// *** Filter condition input format.
	SupBubInfo << "Print out the superbubbles' information:\n";
	SupBubInfo << "Total number of superbubbles:" << setnum << endl;
//	ofstream NewickTree("newicktree.txt");
//	string newicktree = DrawBubTree(affinfo,pairnodeset);
//	NewickTree << newicktree << endl;
	for(int bub = 0;bub < setnum;bub++){
		vector<ssize_t> pairnodeinfo = pairnodeset[bub];
//		if(pairnodeinfo[0] != 1)
//			continue;
		cout << "bubble:" << bub+1 << endl;
		//cout << "pairnodeinfo:" << pairnodeinfo[0] << ";" << pairnodeinfo[1] << ";" << pairnodeinfo[2] << endl;
		vector<int> affiliation = affinfo[bub];
		string bubcolor = bubblecolors[bub];
//		ssize_t start_node = pairnodeinfo[1];
		int level = pairnodeinfo[0];
		BubInfo bubbleinfo;
		vector<int> colorcount(2*numcolors,0);
//		string color = getnodecolor(start_node);
		int affiliationlen = affiliation.size();
		vector<int> visited(dbg.num_nodes(),0);
		bubbleinfo.visited = visited;
		bubbleinfo.level = level;
		bubbleinfo.bubbleid = bub+1;
		bubbleinfo.pairnodeset = pairnodeinfo;
		bubbleinfo.affiliation = affiliation;
		vector<string> branches;
		branches.clear();
		bubbleinfo.branches = branches;
		bubbleinfo.numnodes = 1;
		bubbleinfo.numsupernodes = 1;
		bubbleinfo.numedges = 0;
		bubbleinfo.numbranches = 0;
//		bubbleinfo.firstnode = start_node;
		bubbleinfo.bubcolor = bubcolor;
		bubbleinfo.nodespos = nodepos;
//cout << start_node << endl;
//		BubInfo bubbleinfo_result = get_Bub_Info(nestlevelset,bubbleinfo);
		cout << "start to get c-su detail information" << endl;
		BubInfo bubbleinfo_result = BFSget_Bub_Info(bubbleinfo);//,pairnodelist);
		//BubInfo bubbleinfo_result = BFSget_Bub_Info(nestlevelset,bubbleinfo);//,pairnodelist);
		bubbleinfoset.push_back(bubbleinfo_result);
//		cout << "Superbubble " << groupid << endl;

///############################# Start to output the result ############################

//###### Basic information: bubble id,bubble level,bubble type,father/child bubble id,source/sink node,isolated/contiguous information,nodes,edges,supernodes,branches ######
		SupBubInfo << "Superbubble " << bub+1 << endl;
		SupBubInfo << "  level:" << level << endl;
		SupBubInfo << "  bubblecolor:" << bubcolor << endl;
//		int affiliationlen = affiliation.size();
//		bool simplebubble = false;
		
//		if(affiliationlen == 3 && affiliation[2] == -1 && length == 3 && getnodecolor(pairnodeinfo[1]) == getnodecolor(pairnodeinfo[2]))
//			SupBubInfo << "  Bubble type:Simple bubble" << endl;
//		else{
//			SupBubInfo << "  Bubble type:Complex bubble" << endl;
//		}
//			simplebubble = true;
		if(affiliation[0] != bub+1){
			cout << bub+1 << ":" << affiliation[0] << endl;
			cout << "Oh,my god! There are errors in affiliation's bubble id\n";
		}
		else{
			SupBubInfo << "  father bubble id:" << affiliation[1] << endl;
			SupBubInfo << "  child bubble id:" << affiliation[2];
			for(int k = 3;k<affiliationlen;k++){
				SupBubInfo << ";" << affiliation[k];
			}
			SupBubInfo << endl;
		}
//		int bubchildnum = 0;
//		int bubmaxlevel = pairnodeinfo[0];
//		if(affiliation[2] != -1){
//			bubchildnum = affiliation.size()-2;
//			for(int i = 0;i< bubchildnum;i++){
//				int childid = affiliation[i+2]-1;
//				vector<ssize_t> tmppairnodeinfo = pairnodeset[childid];
//				if(tmppairnodeinfo[0] > bubmaxlevel){
//					bubmaxlevel = tmppairnodeinfo[0];
//				}
//			}
//		}
//		SupBubInfo << "  bubchild number:" << bubchildnum << endl;
//		SupBubInfo << "  bubmaxlevel:" << bubmaxlevel << endl;
		string color1 = getnodecolor(pairnodeinfo[1]);
		color1 = ColorCap({color1,researched_samples});
		string color2 = getnodecolor(pairnodeinfo[2]);
		color2 = ColorCap({color2,researched_samples});
		SupBubInfo << "  source node:" << pairnodeinfo[1] << "\tposition:" << nodepos[pairnodeinfo[1]] << "\tseq:" << dbg.node_label(pairnodeinfo[1]) << "\n\tcolor:" << color1 << endl;
		SupBubInfo << "  sink node:  " << pairnodeinfo[2] << "\tposition:" << nodepos[pairnodeinfo[2]] << "\tseq:" << dbg.node_label(pairnodeinfo[2]) << "\n\tcolor:" << color2 << endl;
		SupBubInfo << "  node number:" << bubbleinfo_result.numnodes << endl;
		SupBubInfo << "  edge number:" << bubbleinfo_result.numedges << endl;
		SupBubInfo << "  supernode number:" << bubbleinfo_result.numsupernodes << endl;
		SupBubInfo << "  branch number:" << bubbleinfo_result.numbranches << endl;
		double double_numnodes = bubbleinfo_result.numnodes;
		double density = (bubbleinfo_result.numedges-double_numnodes+1)/(3*double_numnodes+1);
		SupBubInfo << "  density:" << density << endl;
//		int samplelength = bubbleinfo_result.samplelength;
//		cout << "actual length: " << samplelength << endl;
		int maxlength = bubbleinfo_result.maxlength;
		SupBubInfo << "  longest sequence length: " << maxlength << endl;
//		vector<int> gapbasisid = bubbleinfo_result.gapbasisid;
//		if(gapbasisid.size() > 0)
//			SupBubInfo << "  gap basis id:" << gapbasisid << endl;
/*			
		string variantspos = bubbleinfo_result.variantspos;
		SupBubInfo << "  variant position:" << variantspos << endl;
		int varnum = 0;
		int varlen = variantspos.size();
		for(int j = 0;j < varlen;j++){
			if(variantspos[j] != '0')
		//	if(variantspos[j] == '1')
				varnum++;
		}
		SupBubInfo << "  variant possible number:" << varnum << endl;
*/		
		vector<string> pathseqs = bubbleinfo_result.pathseqs;
		vector<string> pathcolors = bubbleinfo_result.pathcolors;
		vector<int> pathweights = bubbleinfo_result.pathweights;
		vector<int> pathlengths = bubbleinfo_result.pathlengths;
		vector<int> pathsupernodes = bubbleinfo_result.pathsupernodes;
		int pathnum = pathseqs.size();
		for(int i = 0;i < pathnum;i++){
			int pathcolor = 0;
			string tmpstr = pathcolors[i];
			for(int j = 0;j < numcolors;j++){
				if(tmpstr[j] == '1')
					pathcolor++;
			}
			//bubbleid;pathid;pathwidth;pathlength;pathsupernodes;pathweights;path sequence;pathcolors
			Pathinfo << bub+1 << "\t" << i << "\t" << pathcolor << "\t" << pathlengths[i] << "\t" << pathsupernodes[i] << "\t" << pathcolor*pathlengths[i] << "\t" << pathweights[i] << "\t" << pathseqs[i] << "\t" << pathcolors[i] << endl;
			
		}

		vector<vector<int>> pathbranchwidth = bubbleinfo_result.pathbranchwidths;
		Branchwidth << "bubble:" << bub+1 << endl;
		for(int i = 0;i < pathnum;i++){
			Branchwidth << pathbranchwidth[i] << endl;
			
		}
		

		vector<vector<int>> weightmat = bubbleinfo_result.weightmat;
		weightinfo << "bubble:" << bub+1 << "\tbasis length:" << weightmat[0].size() << endl;
		for(int j = 0;j < numcolors;j++){
			weightinfo << weightmat[j] << endl;
		}

		vector<Basis> bubblebasis = bubbleinfo_result.bubblebasis;
		int basisnum = bubblebasis.size();
		cout << "bubble " << bub + 1<< " has " << basisnum << " basis" << endl;
		int bubbasisid = 0;
		for(int j = 0;j < basisnum;j++){
			Basis basis = bubblebasis[j];
			string incomingnodestr = to_string(basis.incomingnode[0]);
			string outgoingnodestr = to_string(basis.outgoingnode[0]);
			if(basis.type == 1){
				int incomingnodenum = basis.incomingnode.size();
				int outgoingnodenum = basis.outgoingnode.size();
				for(int j1 = 1;j1 < incomingnodenum;j1++){
					string tmpincomingnodestr = to_string(basis.incomingnode[j1]);
					incomingnodestr = incomingnodestr+';'+tmpincomingnodestr;
				}
				for(int j1 = 1;j1 < outgoingnodenum;j1++){
					string tmpoutgoingnodestr = to_string(basis.outgoingnode[j1]);
					outgoingnodestr = outgoingnodestr+';'+tmpoutgoingnodestr;
				}
			}
			basisinfo << bub+1 << "\t" << bubbasisid << "\t" << basis.type << "\t" << basis.width << "\t" << basis.length << "\t" << basis.startnode << "\t" << basis.endnode << "\t";
			basisinfo << incomingnodestr << "\t" << outgoingnodestr << "\t" << basis.sequence << endl;
			bubbasisid++;
		}
		vector<int> sampleweight = bubbleinfo_result.sampleweight;
//		for(int j = 0;j < numcolors;j++){
			sampleweightinfo << sampleweight << endl;
		vector<int> samplelength = bubbleinfo_result.samplelength;
		SupBubInfo << "  Samples path length:";
		SupBubInfo << samplelength << endl;

//####### Branch information #######
		vector<string> branchset = bubbleinfo_result.branches;
		vector<vector<ssize_t>> branchnodes = bubbleinfo_result.branchnodes;
		vector<vector<string>> branchcolors = bubbleinfo_result.branchcolors;
		int len = branchset.size();
		string colorstr;
		if(len == 2 && bubcolor == researched_samples){
			for(int i = 0;i < len;i++){
				vector<string> branchcolor = branchcolors[i];
				int branchlen = branchcolor.size();
				string colstr;
				if(branchlen > 2){
					colstr = branchcolor[1];
				}
				else{
					cout << bub+1 << " meet special SNP condition" << endl;
					colstr = branchcolor[1];
				}
				for(int j = 0;j < numcolors;j++){
					vector<int> tmpvec = snpvec[j];
					if(colstr[j] == '1'){
						tmpvec.push_back(1);
					}
					else{tmpvec.push_back(0);}
					snpvec[j] = tmpvec;
				}
			}
		}
		for(int i = 0;i < len;i++){
			vector<ssize_t> branchnode = branchnodes[i];
			int branchsize = branchnode.size();
			SupBubInfo << "  branch" << i+1 << "\n    startnode:" << branchnode[0] << "\tseq:" << dbg.node_label(branchnode[0]) << endl;
			SupBubInfo << "    endnode:" << branchnode[branchsize-1] << "\tseq:" << dbg.node_label(branchnode[branchsize-1]) << endl;
			SupBubInfo << "    length:" << branchsize << endl;
///print out the branch's all possible colors.
			SupBubInfo << "    branchcolor:" << endl;
			vector<string> branchcolor = branchcolors[i];
			colorstr = branchcolor[0];
			int colorsize = branchcolor.size();
		//	cout << branchsize << "\t" << colorsize << endl;
			int colorchangepos = 0;
			for(int j = 1;j < colorsize;j++){
				if(branchcolor[j-1] != branchcolor[j]){
					int num = j-1-colorchangepos+1;
					if(num == 1){
						SupBubInfo << "\tnode:" << branchnode[colorchangepos] << "\tcolor:" << branchcolor[j-1] << "\tnumber:" << num << "\n";
					}
					else{
						SupBubInfo << "\tnode:" << branchnode[colorchangepos] << "->" << branchnode[j-1] << "\tcolor:" << branchcolor[j-1] << "\tnumber:" << num << "\n";
					}
					colorchangepos = j;
				}
				if(j == colorsize-1){
					int num = j-colorchangepos+1;
					if(num == 1){
						SupBubInfo << "\tnode:" << branchnode[colorchangepos] << "\tcolor:" << branchcolor[j] << "\tnumber:" << num << "\n";
					}
					else{
						SupBubInfo << "\tnode:" << branchnode[colorchangepos] << "->" << branchnode[j] << "\tcolor:" << branchcolor[j] << "\tnumber:" << num << "\n";
					}

				}
			}
		//	branchcolor.erase( unique( branchcolor.begin(), branchcolor.end() ), branchcolor.end() );
		//	int uniqcolor = branchcolor.size();
			SupBubInfo << "  sequence:" << branchset[i] << endl;
		}		
		SupBubInfo << "  getend" << endl;
	}
	for(int i = 0;i < numcolors;i++){
		simplesnpinfo << snpvec[i] << endl;
	}
	cout << "Superbubbles' details are listed in the superbubbleinfo.txt\n";
	int totaltime = getMilliSpan(t);
	if(totaltime < 1000){
		cout << "Find bubble information's time: " << getMilliSpan(t) << "ms" << endl;
	}
	else{
		int tolsecond = totaltime/1000;
		int hour = tolsecond/3600;
		int second1 = tolsecond%3600;
		int minute = second1/60;
		int second = second1%60;
		cout << "Find bubble information's time: " << hour << "h " << minute << "m " << second << "s " << endl;
	}
	return bubbleinfoset;
}


//To get the neighbour information of node in dbg
void nearby(){ //get the neighbour information of node
	int t = getMilliCount();
	ofstream nearbyfile("nearbyinfo.txt");
	cout << "Start to get the nodes' nearby information ..." << endl;
//for each candidate node in the graph
	ssize_t num_nodes = dbg.num_nodes();
	
	for(ssize_t node = 0;node< num_nodes;node++){
		int indegree=dbg.indegree(node);
		int outdegree=dbg.outdegree(node);
		nearbyfile << "node:" << node << "\tseq:" << dbg.node_label(node)<< "\tcolor:" << getnodecolor(node) << "\n";
//incoming nodes
		if(indegree==0)
			nearbyfile << "the node doesn't have indegree edge\n";
		else{
			nearbyfile << "the indegree of node is:\t"<<indegree<<"\n";
			for (unsigned long x = 0; x < dbg.sigma + 1; x++) { 
				ssize_t edge = dbg.incoming_edge(node, x);
				ssize_t node1 = dbg.incoming(node, x);
				if (edge == -1)
					continue;
				if (node1 == -1)
					continue;
				nearbyfile << "\tfront node ID:" << node1 << "\tseq:" << dbg.node_label(node1)<<"\n";
				nearbyfile << "\tfront edge ID:" << edge << "\tseq:" << dbg.edge_label(edge)<<"\n";
    	    }
		}
//outgoing nodes
    	if(outdegree==0)
    	    nearbyfile << "the node doesn't have outdegree edge\n";
    	else{
    	    nearbyfile << "the outdegree of node is:\t"<< outdegree<<"\n";
			for (unsigned long x = 1; x < dbg.sigma + 1; x++) { 
				ssize_t edge = dbg.outgoing_edge(node, x);
				if (edge == -1)
					continue;
				ssize_t pos = dbg._edge_to_node(edge);
				for(unsigned long x2 = 0;x2 < 5;x2++){
					ssize_t frontnode = dbg.incoming(pos,x2);
					if(frontnode == -1)
						continue;
					if(frontnode == node){
						edge = dbg.incoming_edge(pos,x2);
						break;
					}
				}
				nearbyfile << "\tnext node ID:" << pos << "\tseq:" << dbg.node_label(pos)<<"\n";
				nearbyfile << "\tnext edge ID:" << edge << "\tseq:" << dbg.edge_label(edge)<<"\n";
    	    }
    	}
	}
	cout << "Nodes' information is listed in the nearbyinfo.txt" << endl;
	//int t = getMilliCount();
	int totaltime = getMilliSpan(t);
	if(totaltime < 1000){
		cout << "Find nearby information's time: " << getMilliSpan(t) << "ms" << endl;
	}
	else{
		int tolsecond = totaltime/1000;
		int hour = tolsecond/3600;
		int second1 = tolsecond%3600;
		int minute = second1/60;
		int second = second1%60;
		cout << "Find nearby information's time: " << hour << "h" << minute << "m" << second << "s" << endl;
	}
}


///To get the file showed in Cytoscape
//void Cytoscape(vector<int> nestlevelset,C_Su c_su){//vector<int> myplusnode){
//void Cytoscape(vector<int> myplusnode){
void Cytoscape(){
	ofstream Cytoscapefile("Cytoscape.txt");
	cout << "Start to get the Cytoscape file" << endl;
	Cytoscapefile << "Source node\tTarget node\tedge color\tnodelevel1\tnodelevel2\tnode1color\tnode2color\n";//edge coverage\n";
    sdsl::bit_vector visited = sdsl::bit_vector(dbg.num_nodes(), 0);
//	vector<int> nestlevelset(dbg.num_nodes(),1);
	vector<int> myplusnode(dbg.num_nodes(),2);
//	vector<int> myplusnode = colcovplus.myplusnode;
//	vector<int> nestlevelset = GetnestinglevelDFS(dbg);//method 2 DFS
	ssize_t num_nodes = dbg.num_nodes();
	int num_colors = colors.size()/dbg.size();
	string zerostr(num_colors,'0');
	for(ssize_t node = 0;node < num_nodes;node++){
		if(!visited[node] && dbg.outdegree(node) != 0 && myplusnode[node] >= 1){
			for(unsigned long x = 0;x < dbg.sigma+1;x++){
//the process to get the edge color and node color is similar,the latter need to visit all the incoming edge of the node
				ssize_t edge = dbg.incoming_edge(node,x);
				if(edge == -1)
					continue;
				ssize_t backnode = dbg.incoming(node,x);
				string outstring = getedgecolor(edge);
				if(outstring == zerostr)
					outstring = getnodecolor(node);
//omit the middle nodes in the bridge, node order: backnode2 => backnode1 => (bridge edge) => nextnode2 => nextnode1
				if(dbg.outdegree(node) == 1 && dbg.indegree(node) == 1){
					ssize_t nextnode1 = node;
					ssize_t nextnode2 = node;
					string outstring1 = "0";
					if(dbg.outdegree(node) == 1){
						while(dbg.indegree(nextnode1) == 1 && dbg.outdegree(nextnode1) == 1){
							for(unsigned long x1 = 1;x1 < dbg.sigma + 1;x1++){
								ssize_t nextedge = dbg.outgoing_edge(nextnode1,x1);
								if(nextedge == -1)
									continue;
								//outstring1 = getedgecolor(nextedge);
								outstring1 = getnodecolor(nextnode1);
								visited[nextnode1] = 1;
								nextnode2 = nextnode1;
								nextnode1 = dbg.outgoing(nextnode1,x1);
								break;
							}
						}
//						if(nextnode2 == 75726)
//							cerr << nextnode2 << "->" << nextnode1 << "\t" << outstring1 << endl;
						Cytoscapefile << nextnode2 << "\t" << nextnode1 << "\t" << outstring1 << "\t" << getnodecolor(nextnode2) << "\t" << getnodecolor(nextnode1) << "\n";
					}
					ssize_t backnode1 = node;
					ssize_t backnode2 = backnode;
					string outstring2 = outstring;
					if(dbg.indegree(node) == 1){
						while(dbg.indegree(backnode2) == 1 && dbg.outdegree(backnode2) == 1){
							for(unsigned long x1 = 0;x1 < dbg.sigma + 1;x1++){
								ssize_t foredge = dbg.incoming_edge(backnode2,x1);
								if(foredge == -1)
									continue;
								outstring2 = getedgecolor(foredge);
								if(outstring2 == zerostr)
									outstring2 = getnodecolor(backnode1);
								visited[backnode2] = 1;
								backnode1 = backnode2;
								backnode2 = dbg.incoming(backnode2,x1);
								break;
							}
						}
//						if(backnode2 == 75726)
//							cerr << backnode2 << "->" << backnode1 << "\t" << outstring2 << endl;
						Cytoscapefile << backnode2 << "\t" << backnode1 << "\t" << outstring2 << "\t" << getnodecolor(backnode2) << "\t" << getnodecolor(backnode1) << "\n";
					}
					//outstring1 = getnodecolor(nextnode2);
					Cytoscapefile << backnode1 << "\t" << nextnode2 << "\t" << outstring << "\t" << getnodecolor(backnode1) << "\t" << getnodecolor(nextnode2) << "\n";
					//node order:backnode2=>backnode1=>(bridge edge)=>nextnode2=>nextnode1
				}
				else{
					Cytoscapefile << backnode << "\t" << node << "\t" << outstring << "\t" << getnodecolor(backnode) << "\t" << getnodecolor(node) << "\n";
				}
			}
			visited[node] = 1;
		}
	}
	cout << "Cytoscape file is listed in the Cytoscape.txt" << endl;
}
//############################################################################################


void visited_dbg(){
//	for(int i = 0;i<1024;i++){
//	cout << "sample id:" << i << endl;
	ssize_t startnode = 84931;
	int numnodes = dbg.num_nodes();
	int numcolors = colors.size()/dbg.size();
	string stdcolor(numcolors,'1');
	string refcolor(numcolors-1,'0');
	vector<int> visited(numnodes,0);
	queue<ssize_t> nodelist;
	nodelist.push(startnode);
	int nodecount = 1;
	visited[startnode] = 1;
	while(!nodelist.empty()){
		ssize_t node = nodelist.front();
		for(unsigned long x = 1;x < 5; x++){
			ssize_t nextnode1 = dbg.outgoing(node,x);
			if(nextnode1 == -1 || visited[nextnode1] == 1)
				continue;
			string nodecolor = getnodecolor(nextnode1);
//			if(nodecolor[i]== '0')
//				continue;
			cout << dbg.node_label(nextnode1) << "\t" << node << "->" << nextnode1 << "\t" << nodecolor << endl;
			if(nodecolor == stdcolor){
				cerr << "nextnode:" << nextnode1 << endl;
				break;
			}
			nodelist.push(nextnode1);
			nodecount++;
			visited[nextnode1] = 1;
		}
		nodelist.pop();
	}
	cout << "number of connected nodes: " << nodecount << endl;
//	}
}

string read_researched_samples_from_file(string researched_samples_file){
	string myresearched_samples;
	ifstream myresearchedsamp_file;
	myresearchedsamp_file.open(researched_samples_file.data());
	//string myresearched_sample;
	if(!myresearchedsamp_file.is_open()){
		cout << "open researched file failure!" << endl;
	}
	else{
		while(myresearchedsamp_file >> myresearched_samples){
//			return myresearched_samples;
		}
	}
	int numcolors = colors.size()/dbg.size();
	int len = myresearched_samples.size();
	if(len != numcolors){
		cout << "researched file has problem" << endl;
	}
	return myresearched_samples;
}


/////To get the path of one node around and each node's indegree and outdegree are both 1
void GetPathofOneCutnode(ssize_t cutnode){
    unsigned long count = 1;
    cout << "nearby node of " << cutnode << "\n";
    if(dbg.outdegree(cutnode) == 1){
        for(unsigned long x = 1;x < dbg.sigma + 1;x++){
            ssize_t edge = dbg.outgoing_edge(cutnode,x);
            if(edge == -1)
                continue;
            ssize_t pos = dbg.outgoing(cutnode,x);
            cout << cutnode << "=>" << pos;
            count++;
            while(dbg.outdegree(pos) == 1 && dbg.indegree(pos) == 1){
                ssize_t next_edge = 0;
                for(unsigned long x2 = 1;x2 < dbg.sigma + 1;x2++){
                    next_edge = dbg.outgoing_edge(pos,x2);
                    if(next_edge != -1){
                        pos = dbg.outgoing(pos,x2);
                        break;
                    }
                }
                cout << "=>" << pos;
                count++;
            }
        }
    }
    cout << "\n";
    if(dbg.indegree(cutnode) == 1){
        for(unsigned long x = 0;x < dbg.sigma + 1;x++){
            ssize_t edge = dbg.incoming_edge(cutnode,x);
            if(edge == -1)
                continue;
            ssize_t pos = dbg.incoming(cutnode,x);
            if(pos == -1)
                continue;
            cout << cutnode << "<=" << pos;
            count++;
            while(dbg.indegree(pos) == 1 && dbg.outdegree(pos) == 1){
                ssize_t next_edge = 0;
                for(unsigned long x2 = 0;x2 < dbg.sigma + 1;x2++){
                    next_edge = dbg.incoming_edge(pos,x2);
                    if(next_edge != -1){
                    //  cout << next_edge <<"\n" << dbg.edge_label(next_edge) << "\n";
                        pos = dbg.incoming(pos,x2);
                        break;
                    }
                }
                cout << "<=" << pos;
                count++;
            }
        }
    }
    cout << "\nlength of branch is:" << count << "\n";
    return;
}

///reconstruct the genome sequences
void reconstruct_genome_seq(vector<int> sampleid){
	int numcolors = colors.size()/dbg.size();
	int numsampid = sampleid.size();
	if(numsampid > numcolors){
		cout << "There are too many elements in the sample id set" << endl;
		exit(0);
	}
//	string zerolabel = dbg.node_label(0);
//	int kmer = zerolabel.size();
	int kmer = dbg.k;
	string startlabel = dbgstartstring.substr(0,kmer);
	ssize_t startnode = -1;
	int numnodes = dbg.num_nodes();
	for(ssize_t node = 0;node < numnodes;node++){
		string tmpnodelabel = dbg.node_label(node);
		if(tmpnodelabel == startlabel){
			startnode = node;
		//	cout << "find start node: " << node << endl;
			break;
		}
	}
//	cout << startlabel << endl;
	for(int id = 0;id < numsampid;id++){
		int idvalue = sampleid[id];
		bool getend = false;
		while(getend == false){
			for(int k = 1;k < 5;k++){
				ssize_t nextnode = dbg.outgoing(startnode,k);
				if(nextnode == -1){
					if(k == 4)
						getend = true;
					continue;
				}
				string nodecol = getnodecolor(nextnode);
				if(nodecol[idvalue-1] == '1'){
					startnode = nextnode;
					startlabel += base[k];
				//	cout << startlabel << endl;
					break;
				}
				if(k == 4)
					getend = true;
			}
		}
		cout << ">" << idvalue << "\n" << startlabel << endl;
	}
	return;	
}

void judge_circle(){
	ssize_t numnodes = dbg.num_nodes();
	int numcolors = colors.size()/dbg.size();
	string zerocolor(numcolors,'0');
	vector<string> colorset;
	for(ssize_t node = 0;node < numnodes;node++){
		string circlecol = zerocolor;
		for(unsigned long x2 = 0;x2 < dbg.sigma+1;x2++){
			ssize_t frontnode = dbg.outgoing(node,x2);
			//ssize_t frontnode = dbg.incoming(node,x2);
			if(frontnode == -1)
				continue;

			string tmpnodecolor = getnodecolor(frontnode);
			if(tmpnodecolor != zerocolor)
				colorset.push_back(tmpnodecolor);
		//	int colcmp = StringCMP(circlecol,tmpnodecolor);
		//	if(colcmp == 1 ||(colcmp == 2 && tmpnodecolor != zerocolor)){
		//		cerr << "meet circle at " << node << endl; 
		//	}
		//	else{
		//		circlecol = ColorCup({circlecol,tmpnodecolor});
		//	}
		}
		int num = colorset.size();
		if(num > 1){
			string capcolor = ColorCap(colorset);
			if(capcolor != zerocolor)
				cerr << "meet circle at " << node << endl;
		}
	}
}

void GetAllWeight(ssize_t startnode){
//	ssize_t startnode = 553609;
//	int kmer = dbg.k-1;
	ofstream weightfile("allweight.txt");
    cout << "Starting to get all sample weight directly...\n";
	int numcolors = colors.size()/dbg.size();
	vector<int> visited(dbg.num_nodes(),0);
	visited[startnode] = 1;//0:unvisited 1:visiting 2:visited 3:incoming node unvisited completely
	string stdcolor(numcolors,'1');
	string zerocolor(numcolors,'0');
	vector<int> emptyvec;
	vector<vector<int>> weightmat(numcolors,emptyvec);
	queue<ssize_t> nodelist;
	nodelist.push(startnode);
	int count = 0;
	while(!nodelist.empty()){
		ssize_t supernode = nodelist.front();
		count++;
		string nodecolor = getnodecolor(supernode);
		nodecolor = ColorCap({nodecolor,researched_samples});
		int weight;
		int basiswidth = 0;
		for(int i = 0;i < numcolors;i++){
			if(nodecolor[i] == '1'){
				basiswidth++;
			}
		}
		if(basiswidth == 0)
			cout << supernode << " 's color is zerocolor,error!!"<< endl;
//		weight = kmer * basiswidth;	
		weight = basiswidth;	
		for(int i = 0;i < numcolors;i++){
			vector<int> weightvec = weightmat[i];
			if(nodecolor[i] == '1'){
				weightvec.push_back(weight);
			}
			else{weightvec.push_back(0);}
			weightmat[i] = weightvec;
		}
		for(int x = 1;x <5;x++){
			ssize_t nextnode = dbg.outgoing(supernode,x);
			if(nextnode == -1|| visited[nextnode] == 2)
				continue;
			int branchlen = 0;
			string branchcol;
			while(dbg.outdegree(nextnode)==1 && dbg.indegree(nextnode) == 1){
				visited[nextnode] = 2;
				branchlen++;
				count++;
				if(branchlen == 1){
					branchcol = getnodecolor(nextnode);
					branchcol = ColorCap({branchcol,researched_samples});
				}
				for(int x2 = 1;x2 < 5;x2++){
					ssize_t pos = dbg.outgoing(nextnode,x2);
					if(pos == -1)
						continue;
					if(visited[pos]==2){
						//cout << "visiting meet error"<< endl;
						continue;
					}
					nextnode = pos;
					break;
				}
			}
			if(branchlen > 0){
				vector<int> coloredsamp1;
				for(int i = 0;i < numcolors;i++){
					if(branchcol[i] == '1')
						coloredsamp1.push_back(i);
				}
				weight = branchlen * coloredsamp1.size();	
				for(int i = 0;i < numcolors;i++){
					vector<int> weightvec = weightmat[i];
					if(branchcol[i] == '1'){
						weightvec.push_back(weight);
					}
					else{weightvec.push_back(0);}
					weightmat[i] = weightvec;
				}
			}
			if(visited[nextnode] == 0){
				nodelist.push(nextnode);
				visited[nextnode] = 1;
			}
		}
		visited[supernode] = 2;
		nodelist.pop();
		cout << count << endl;
	}
	for(int i = 0;i < numcolors;i++){
		weightfile << weightmat[i] << endl;
	}
}


void GetAllWeight2(){
	ssize_t startnode = 37614;
//	ssize_t startnode = 88965;
//	int kmer = dbg.k-1;
	ofstream weightfile("allweight2.txt");
    cout << "Starting to get all sample weight directly...\n";
	int numcolors = colors.size()/dbg.size();
	vector<int> visited(dbg.num_nodes(),0);
	visited[startnode] = 1;//0:unvisited 1:visiting 2:visited 3:incoming node unvisited completely
	string stdcolor(numcolors,'1');
	string zerocolor(numcolors,'0');
	vector<int> emptyvec;
	vector<vector<int>> weightmat(numcolors,emptyvec);
	queue<ssize_t> nodelist;
	nodelist.push(startnode);
	int nodecount = 0;
	while(!nodelist.empty()){
		ssize_t supernode = nodelist.front();
		nodecount++;
		cout << nodecount << endl;
		string nodecolor = getnodecolor(supernode);
		nodecolor = ColorCap({nodecolor,researched_samples});
//		int weight;
//		int basiswidth = 0;
//		for(int i = 0;i < numcolors;i++){
//			if(nodecolor[i] == '1'){
//				basiswidth++;
//			}
//		}
//		if(basiswidth == 0)
//			cout << supernode << " 's color is zerocolor,error!!"<< endl;
//		weight = basiswidth;	
		if(nodecolor != zerocolor && nodecolor != researched_samples){
		for(int i = 0;i < numcolors;i++){
			vector<int> weightvec = weightmat[i];
			if(nodecolor[i] == '1'){
		//		weightvec.push_back(weight);
				weightvec.push_back(1);
			}
			else{weightvec.push_back(0);}
			weightmat[i] = weightvec;
		}
		}
		for(int x = 1;x <5;x++){
			ssize_t nextnode = dbg.outgoing(supernode,x);
			if(nextnode == -1|| visited[nextnode] == 2)
				continue;
			//int branchlen = 0;
			string branchcol;
		//	bool end = false;
			while(dbg.outdegree(nextnode)==1 && dbg.indegree(nextnode) == 1){
				nodecount++;
				cout << nodecount << "\t" << nextnode << endl;
				visited[nextnode] = 2;
			//	branchlen++;
				branchcol = getnodecolor(nextnode);
				branchcol = ColorCap({branchcol,researched_samples});
			//	vector<int> coloredsamp1;
	//			int width = 0;
	//			for(int i = 0;i < numcolors;i++){
	//				if(branchcol[i] == '1')
	//					width++;
	//			}
	//			weight = width;	
		if(nodecolor != zerocolor && nodecolor != researched_samples){
				for(int i = 0;i < numcolors;i++){
					vector<int> weightvec = weightmat[i];
					if(branchcol[i] == '1'){
				//		weightvec.push_back(weight);
						weightvec.push_back(1);
					}
					else{weightvec.push_back(0);}
					weightmat[i] = weightvec;
				}
		}
			//	if(branchlen == 1)
			//		branchcol = getnodecolor(nextnode);
				for(int x2 = 1;x2 < 5;x2++){
					ssize_t pos = dbg.outgoing(nextnode,x2);
					if(pos == -1)
						continue;
			//		if(visited[pos]==2){
			//			//cout << "visiting meet error"<< endl;
			//			if(x2 == 4)
			//				end = true;
			//			else{
			//				continue;
			//			}
			//		}
					nextnode = pos;
					break;
				}
		//		if(end == true)
		//			break;
//				cout << nextnode << endl;
			}
			if(visited[nextnode] == 0){
				nodelist.push(nextnode);
				visited[nextnode] = 1;
			}
		}
		visited[supernode] = 2;
		nodelist.pop();
	}
	for(int i = 0;i < numcolors;i++){
		weightfile << weightmat[i] << endl;
	}
}


void visiting_samp_nodes(vector<int> samps,ssize_t startnode){
	int samp_num = samps.size();
	for(int i = 0;i < samp_num;i++){
		int selectid = samps[i];
		cout << "selected sample id:" << selectid << endl;
//		bool getend = false;
//		bool getstart = false;
		cout << startnode;
		int nodenum = 1;
		ssize_t pos = startnode;
		while(1){
			int outgoingnum = 0;
			ssize_t tmpnode = -1;
			for(int x = 1;x < 5;x++){
				ssize_t nextnode = dbg.outgoing(pos,x);
				if(nextnode == -1)
					continue;
				string nodecolor = getnodecolor(nextnode);
				if(nodecolor[selectid - 1] == '1'){
					cout << "->" << nextnode;
					tmpnode = nextnode;
					outgoingnum++;
					nodenum++;
				}
			}
			if(outgoingnum > 1){
				cout << "\nnode " << pos << " has two exports\n";
				break;

			}
			if(tmpnode == -1){
				cout << "\nnode " << pos << " has gotten end\n";
				break;
			}
			pos = tmpnode;
		}
		cout << startnode;
		pos = startnode;
		while(1){
		//	ssize_t pos = startnode;
			int incomingnum = 0;
			ssize_t tmpnode = -1;
			for(int x = 1;x < 5;x++){
				ssize_t frontnode = dbg.incoming(pos,x);
				if(frontnode == -1)
					continue;
				string nodecolor = getnodecolor(frontnode);
				if(nodecolor[selectid - 1] == '1'){
					cout << "<-" << frontnode;
					tmpnode = frontnode;
					incomingnum++;
					nodenum++;
				}
			}
			if(incomingnum > 1){
				cout << "\nnode " << pos << " has two entrances\n";
				break;

			}
			if(tmpnode == -1){
				cout << "\nnode " << pos << " has gotten start position\n";
				break;
			}
			pos = tmpnode;
		}
		cout << "find " << nodenum << " nodes in the path" << endl;
	}
}


//############################################################################################
int main(int argc, char* argv[]) {
//ld -stack=0x40000000;
//	system("ulimit -s 1024000000");	
	int t = getMilliCount();
  parameters_t p;
  parse_arguments(argc, argv, p);
  cout << "pack-color compiled with supported colors=" << NUM_COLS << std::endl;
  //ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
  // Can add this to save a couple seconds off traversal - not really worth it.
  cout << "loading dbg" << std::endl;
  //debruijn_graph_shifted<> dbg;
  load_from_file(dbg, p.input_filename);
  //input.close();
 // cout << "loading colors" << std::endl;
 // sd_vector<> colors;
  load_from_file(colors, p.color_filename);

  cout << "k             : " << dbg.k << endl;
  cout << "num_nodes()   : " << dbg.num_nodes() << endl;
  cout << "num_edges()   : " << dbg.num_edges() << endl;
  cout << "colors        : " << colors.size() / dbg.size() << endl; 
  cout << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
  cout << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;
  cout << "Color size    : " << size_in_mega_bytes(colors) << " MB" << endl;

//	vector<int> plusnodes = GetPlusnode();
//	vector<int> plusnodes(dbg.num_nodes(),1);
	int numcolors = colors.size()/dbg.size();
	string zerocolor(numcolors,'0');
	researched_samples.assign(numcolors,'1');
//	string researched_samples_file = "researched_samples.txt";
//	researched_samples = read_researched_samples_from_file(researched_samples_file);
	cout << researched_samples << endl;
//	for(int i =0;i <10;i++){
//		researched_samples[i] = '0';
//	}
//	visiting_samp_nodes({8},13163);
//	nearby();
//	C_Su c_su = GetColoredSuperbubble();

	C_Su c_su = read_mypairs_from_file("mypairnodes.txt");
//	vector<vector<ssize_t>> pairnodeset = c_su.pairnodesets;
//	vector<string> bubblecolor = c_su.bubblecolors;
//	vector<ssize_t> pairnodes = pairnodeset[1];
//	cout << "pairnodes:" << pairnodes[0] << "\t" << pairnodes[1] << ":" << pairnodes[2] << "\t" << bubblecolor[1] << endl;
	
//	vector<BubInfo> bubbleinfoset = SupperBubblesAnalysis(c_su);
//	Cytoscape();
//	GetPathofOneCutnode(91879);
//	visited_dbg();
//GetAllWeight(7242376);
//	GetAllWeight2();

	int totaltime = getMilliSpan(t);
	if(totaltime < 1000){
		cout << "Total calculation time: " << getMilliSpan(t) << "ms" << endl;
	}
	else{
		int tolsecond = totaltime/1000;
		int hour = tolsecond/3600;
		int second1 = tolsecond%3600;
		int minute = second1/60;
		int second = second1%60;
		cout << "Total calculation time: " << hour << "h" << minute << "m" << second << "s" << endl;
	}
}
