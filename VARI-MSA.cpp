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
//#include "pg.hpp"

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

static char base[] = {'$','A','C','G','T'};

static debruijn_graph_shifted<> dbg;
static sd_vector<> colors;
static vector<int> mycolors;
//static int REF_LEN = 16569;
static string researched_samples;

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

struct C_Su{
	vector<vector<ssize_t>> pairnodesets;
	vector<vector<int>> affiliationsets;
	vector<string> bubblecolors;
	map<ssize_t,int> nodepos;
	map<ssize_t,int> refnodepos;
	map<int,int> bub2gap;
	map<int,int> bub2var;
};


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

//####################################################################################################################################################
//########################################################  Multiple Sequence Alignment ##############################################################
//####################################################################################################################################################
struct Seqalign_result{
	string alignseq;
	vector<int> basescore;
};

//Seqalign_result SeqAlign(string str,string ref,vector<int> basescore,map<int,char> tmppos2keybase){
Seqalign_result SeqAlign(string str,string ref,vector<int> basescore,map<int,char> tmppos2keybase){
	int m = str.size();
	int n = ref.size();
	int score[(m+1)*(n-m+1)];
	int source[(m+1)*(n-m+1)];
	int range = n - m + 1;
	int keybasevalue = 2;
	Seqalign_result seqalign;
//	vector<int> basescore(n,1);
	score[0] = 0;
	source[0] = -1;
	for(int i = 1;i < n-m+1;i++){
		score[i] = 0;
		source[i] = 1;
	}
	for(int i = 1;i < m+1;i++){
//		cout << "i:" << i << " " << str[i-1] << " " << ref[i-1] << " " << tmppos2keybase[i-1] << endl;
//		seq2samples.find(tmpstr) != seq2samples.end()
		if(str[i-1] == ref[i-1]){
			bool isacgt = false;
			if(tmppos2keybase[i-1] == 'A'||tmppos2keybase[i-1] == 'C'||tmppos2keybase[i-1] == 'G'||tmppos2keybase[i-1] == 'T')
				isacgt = true;
			if(isacgt == true && tmppos2keybase[i-1] == str[i-1]){
//			if(tmppos2keybase.find(i-1) != tmppos2keybase.end() && tmppos2keybase[i-1] == str[i-1])
				score[i*range] = score[(i-1)*range]+keybasevalue;
			}else{
//				cout << i-1 << ";" << basescore[i-1] << endl;
				score[i*range] = score[(i-1)*range]+basescore[i-1];
			}
//			score[i*range] = score[(i-1)*range]+basescore[i-1];
//			basescore[i-1]++;
			source[i*range] = 3;
		}
		else{
			score[i*range] = score[(i-1)*range];
			source[i*range] = 0;
//			if(basescore[i-1] > 1)
//				basescore[i-1]--;
		}
//		cout << i << "\t" << source[i*range] << "\t" << score[i*range]  << endl;
	}
	for(int i = 1;i < m + 1;i++){
		for(int j = 1;j < n - m + 1;j++){
			int tmpscore = score[i*range+j-1];
			int score1 = score[(i-1)*range+j];
			bool diag = false;
			if(str[i-1] == ref[i + j - 1]){
				bool isacgt = false;
				if(tmppos2keybase[i+j-1] == 'A'||tmppos2keybase[i+j-1] == 'C'||tmppos2keybase[i+j-1] == 'G'||tmppos2keybase[i+j-1] == 'T')
					isacgt = true;
			//	if(isacgt == true && tmppos2keybase[i-1] == str[i-1])
				if(isacgt == true && tmppos2keybase[i+j-1] == str[i-1]){
					score1 = score[(i-1)*range+j]+keybasevalue;
				}else{
					score1 = score[(i-1)*range+j]+basescore[i-1];
				}
				diag = true;
			}
			if(tmpscore > score1){
				score[i*range+j] = tmpscore;
				source[i*range+j] = 1;
//				cout << i << ":" << j << "\t" << source[i*range+j] << endl;
			}
			else if(tmpscore < score1){
				score[i*range+j] = score1;
				source[i*range+j] = 0;
				if(diag == true){
					source[i*range+j] = 3;
				}
//				cout << i << ":" << j << "\t" << source[i*range+j] << endl;
			}
			else{
				score[i*range+j] = score1;
				source[i*range+j] = 2;
				if(diag == true){
					source[i*range+j] = 4;
				}
//				cout << i << ":" << j << "\t" << source[i*range+j] << endl;
			}
//			cout << i << ":" << j << "\t" << source[i*range+j] << "\t" << score[i*range+j]  << endl;
		}
	}
//	cout << ref << "\n" << str << endl;
//	for(int i = 0;i < m + 1;i++){
//		for(int j = 0;j < n - m + 1;j++){
//			cout << score[i*range+j] << "(" << source[i*range+j] << ") ";
//
//
//		}
//		cout << endl;
//	}
//	int maxscore = score[m*(n-m)];
	bool getend = false;
	int row = m;
	int col = n-m;
	int pos = n;
	int rowpos = m;
	string alignedseq(n,'-');
//	cout << "range:" << range << endl;
	
	while(getend == false){
		int tmpsource = source[row*range+col];
//		cout << row << "\t" << col << "\t" << tmpsource << "\tscore:" << score[row*range+col] << endl;
		if(tmpsource == 0 || tmpsource == 2||tmpsource == 3||tmpsource == 4){
			row--;
	//		col--;
			alignedseq[pos-1] = str[rowpos-1];
			pos--;
			rowpos--;
//			if(tmpsource == 3){
//				basescore[pos]++;
//			}
//			else{
//				if(basescore[pos] > 1)
//					basescore[pos]--;
//			}
		}
		else if(tmpsource == 1){
			col--;
			pos--;	
//			if(basescore[pos] > 1)
//				basescore[pos]--;
		}
		else{
			getend = true;
		}
	}
	seqalign.alignseq = alignedseq;
	seqalign.basescore = basescore;
//	return alignedseq;
	return seqalign;

}


struct Seqalign1_result{
	vector<string> alignseq;
	map<int,char> pos2keybase;
	vector<int> basescore;
};

Seqalign1_result SeqAlign1(string str,string ref,map<int,char> tmppos2keybase,vector<int> basescore,int maxgap = 2,int multiple = 3){//是否也存在gap，如果有，进行添加
	int n = str.size();
	int reflen = ref.size();
	if(n != reflen){
		Seqalign_result sa = SeqAlign(str,ref,basescore,tmppos2keybase);
		str = sa.alignseq;
		//cout << "align step1:\n" << str << endl;
		basescore = sa.basescore;
	}
	n = str.size();
	int diffcount = 0;
	for(int i = 0;i < n;i++){
		if(ref[i] != str[i])
			diffcount++;
	}
	double rate = double(diffcount)/double(n);
	if(rate > 0.3 || diffcount > 5){
	//range = 3;
		int range = 2*maxgap+1;
		int score[(n+1)*range];
		int source[(n+1)*range];
	//	int multiple = 3;
		int alignscore = 1;
		int snpscore = 1;
		int gapscore = (alignscore + snpscore)*multiple - alignscore - 1;
	//	cout << "gapscore:" << gapscore << endl;
		for(int i1 = 1;i1 < n + 1;i1++){//vertical, initialization
			score[i1*range] = -i1*gapscore;
			source[i1*range] = 2;
		}	
		for(int i2 = 1;i2 < maxgap + 1;i2++){//horizontal
			score[i2] = -i2*gapscore;
			source[i2] = 1;
		}
		score[0] = 0;
		source[0] = -1;
		for(int i1 = 0;i1 < n+1;i1++){
			for(int i2 = 0;i2 < range;i2++){
				if(i1 + i2 > n + maxgap || i2 - i1 > maxgap){//bottom right & top right
					source[i1*range+i2] = -2;
					score[i1*range+i2] = -(gapscore*n);
				}
			}
		}
		for(int i = 1;i < n + 1;i++){
		//	int j0 = i > n - maxgap ? i - n + maxgap : 1;
			for(int j = 0;j < range;j++){
				if(i + j > n + maxgap || j - i > maxgap){continue;}
				if(i < maxgap + 1 && j > 0){//0:diag 1:horizontal 2:vertical 3:1+2
					int score1 = score[(i-1)*range + j-1];
					if(str[i-1] == ref[j-1]){
						bool isacgt = false;
						if(tmppos2keybase[i-1] == 'A'||tmppos2keybase[i-1] == 'C'||tmppos2keybase[i-1] == 'G'||tmppos2keybase[i-1] == 'T')
							isacgt = true;
						if(isacgt == true && tmppos2keybase[i-1] == str[i-1]){
							score1 = score1 + 10;
						}else{
							score1 += alignscore;
						}
					}else{
						score1 = score[(i-1) * range + j-1] - snpscore;
					}
					
					int score2 = score[(i-1)*range+j] - gapscore;// vertical
					int score3 = score[i*range+j-1] - gapscore;// horizontal
				//	cout << i << "\t" << j << "\n" << str[i-1] << "\t" << ref[j-1] << "\n" << score1 << "\t" << score2 << "\t" << score3 << endl;
					if(score1 > score2){
						if(score1 >= score3){
							score[i*range+j] = score1;
							source[i*range+j] = 0;
						}
						else{
							score[i*range+j] = score3;
							source[i*range+j] = 1;
						}
					}
					else if(score1 < score2){
						if(score2 > score3){
							score[i*range+j] = score2;
							source[i*range+j] = 2;
						}
						else if(score2 < score3){
							score[i*range+j] = score3;
							source[i*range + j] = 1;
						}
						else{
							score[i*range +j] = score3;
							source[i*range + j] = 3;
						}
					}
					else{
						if(score1 >= score3){
							score[i*range+j] = score1;
							source[i*range + j] = 0;// diag priority
						}
						else{
							score[i*range+j] = score3;
							source[i*range+j] = 1;
						}
					}
		//			cout << source[i*range+j] << "\t" << score[i*range+j] << endl;
				}
				else if(i >= maxgap + 1){//1:horizontal 4:vertical(equal to 0 before) 5:topright direction(equal to 2 before) 6:1+5
					int raw_j = i + j - maxgap;
		//			cout << "raw_j:"  << i << ":" << j << "\t" << raw_j << endl;
					if(j == range - 1){
						//int score1 = score[(i-1)*range + j] + alignscore;
						int score1 = score[(i-1)*range + j];
						if(str[i-1] == ref[raw_j-1]){
							bool isacgt = false;
							if(tmppos2keybase[i-1] == 'A'||tmppos2keybase[i-1] == 'C'||tmppos2keybase[i-1] == 'G'||tmppos2keybase[i-1] == 'T')
								isacgt = true;
							if(isacgt == true && tmppos2keybase[i-1] == str[i-1]){
								score1 = score1 + 10;
							}else{
								score1 += alignscore;
							}
						}else{
					//		if(str[i-1] != ref[raw_j-1]){
							score1 = score[(i-1) * range + j] - snpscore;
						}				
						int score3 = score[i*range+j-1] - gapscore;// horizontal
						//int score2 = score[(i-1)*range+j] - gapscore;// vertical
						if(score1 >= score3){
							score[i*range+j] = score1;
							source[i*range+j] = 4;
						}
						else{
							score[i*range+j] = score3;
							source[i*range+j] = 1;
						}
					}
					else{//
						//int score1 = score[(i-1)*range + j] + alignscore;
						int score1 = score[(i-1)*range + j];
						if(str[i-1] == ref[raw_j-1]){
							bool isacgt = false;
							if(tmppos2keybase[i-1] == 'A'||tmppos2keybase[i-1] == 'C'||tmppos2keybase[i-1] == 'G'||tmppos2keybase[i-1] == 'T')
								isacgt = true;
							if(isacgt == true && tmppos2keybase[i-1] == str[i-1]){
								score1 = score1 + 10;
							}else{
								score1 += alignscore;
							}
						}else{
					//		if(str[i-1] != ref[raw_j - 1]){
							score1 = score[(i-1) * range + j] - snpscore;
						}				
						int score2 = score[(i-1)*range+j+1] - gapscore;// vertical
						int score3 = score[i*range+j-1] - gapscore;// horizontal
		//			cout << i << "\t" << j << "\n" << str[i-1] << "\t" << ref[raw_j-1] << "\n" << score1 << "\t" << score2 << "\t" << score3 << endl;
						if(score1 > score2){
							if(score1 >= score3){
								score[i*range+j] = score1;
								source[i*range+j] = 4;
							}
							else{
								score[i*range+j] = score3;
								source[i*range+j] = 1;
							}
						}
						else if(score1 < score2){
							if(score2 > score3){
								score[i*range+j] = score2;
								source[i*range+j] = 5;
							}
							else if(score2 < score3){
								score[i*range+j] = score3;
								source[i*range + j] = 1;
							}
							else{
								score[i*range +j] = score3;
								source[i*range + j] = 6;
							}
						}
						else{
							if(score1 >= score3){
								score[i*range+j] = score1;
								source[i*range + j] = 4;// diag priority
							}
							else{
								score[i*range+j] = score3;
								source[i*range+j] = 1;
							}
						}
					}
		//			cout << source[i*range+j] << "\t" << score[i*range+j] << endl;
				}
			}
		}
		bool getend = false;
		string str1 = "";
		string str2 = "";		
		int row = n;
		int col = maxgap;
		int maxlen = n;
	//	int colpos = n;
	//	int rowpos = n;
	//	cout << str << "\t" << ref << endl;
		while(getend == false){
			int tmpsource = source[row*range+col];
	//		cout << "tmpsource:" << tmpsource << endl;
	//		cout << str1 << "\n" << str2 << endl;
	//		cout << row << "\t" << col << endl;
			if(row < maxgap + 1){
				if(tmpsource == 0){
					str1 = str[row-1] + str1;
					str2 = ref[col-1] + str2;
					col--;
					row--;
				}
				else if(tmpsource == 1||tmpsource == 3){
					str1 = "-"+ str1;
					str2 = ref[col-1] + str2;
					col--;
				}
				else if(tmpsource == 2){
					str1 = str[row-1] + str1;
					str2 = "-" + str2;
					for(int i0 = maxlen;i0 >= row;i0--){
						if(tmppos2keybase[i0]){
							tmppos2keybase[i0+1] = tmppos2keybase[i0];
							tmppos2keybase.erase(i0);
						}
					}
					maxlen++;
					row--;
				}
				else{
					getend = true;
				}
			}
			else{
				int raw_col = row + col - maxgap;
				if(tmpsource == 4){
					str1 = str[row-1] + str1;
					str2 = ref[raw_col - 1] + str2;
					row--;
				}
				else if(tmpsource == 1||tmpsource == 6){
					str1 = "-"+ str1;
					str2 = ref[raw_col-1] + str2;
					col--;
				}
				else if(tmpsource == 5){
					str1 = str[row-1] + str1;
					str2 = "-" + str2;
					for(int i0 = maxlen;i0 >= row;i0--){
						if(tmppos2keybase[i0]){
							tmppos2keybase[i0+1] = tmppos2keybase[i0];
							tmppos2keybase.erase(i0);
						}
					}
					maxlen++;
					row--;
					col++;
				}
				else{
					getend = true;
				}
			}
		}
		ref = str2;
		str = str1;
	}
	Seqalign1_result result1;
	//result1.alignseq = {str2,str1};
	result1.alignseq = {ref,str};
	result1.pos2keybase = tmppos2keybase;
	result1.basescore = basescore;
	return result1;
//	return{str2,str1};
}


struct MSA{
	vector<vector<ssize_t>> pairnodesets;
	vector<int> bubid;
	map<ssize_t,int> nodepos;
	vector<string> bubblecolors;
	vector<vector<string>> bubbleseqs;
	vector<vector<int>> bubblelens;
};


MSA GetMSA(C_Su c_su){
	int t = getMilliCount();
	vector<vector<ssize_t>> pairnodeset = c_su.pairnodesets;
	vector<vector<int>> affiliationset = c_su.affiliationsets;
	vector<string> bubblecolors = c_su.bubblecolors;
	map<ssize_t,int> nodepos = c_su.nodepos;
	map<ssize_t,int> refnodepos = c_su.refnodepos;
//	cout << "nodepos:" << nodepos[33745] << "\t" << nodepos[9530] << endl;
	ofstream MyMSA("myMSA.txt");
	ofstream mymergeMSA("mymergeMSA.txt");
	ofstream mymsavar("mymsavar.txt");
	int numcolors = colors.size()/dbg.size();
	string zerocolor(numcolors,'0');
	string stdcolor(numcolors,'1');
	int setnum = pairnodeset.size();
	cout << "bubble number: " << setnum << endl;
	int kmer = dbg.k - 1;
	MSA msa_set;
	msa_set.pairnodesets = pairnodeset;
	msa_set.bubblecolors = bubblecolors;
	msa_set.nodepos = nodepos;
	vector<int> bubid;
	vector<vector<string>> bubbleseqs;
	vector<vector<int>> bubblelens;
	vector<string> msaseqs;
	map<int,char> pos2keybase;
	map<int,double> pos2keybasevalue;
//	string emptystr = "";
	bool firstbub = true;
	bool pre_cutnode = false;
//	int startbub = 0;
	for(int bub = 0;bub < setnum;bub++){
		string bubcolor = bubblecolors[bub];
		if(bubcolor != researched_samples){continue;}
//		string longestsamp = bubcolor;
		vector<ssize_t> pairnodeinfo = pairnodeset[bub];
//		vector<int> aff = affiliationset[bub];
		ssize_t source = pairnodeinfo[1];
		ssize_t sink = pairnodeinfo[2];
		int sourcepos = nodepos[source];
		int sinkpos = nodepos[sink];
		if(firstbub == true){
			string tailbridge = "";
			bool nofathernode = false;
			ssize_t startnode = sink;
			while(nofathernode == false){
				nofathernode = true;
				for(int x = 1;x < 5;x++){
					ssize_t fathnode = dbg.outgoing(startnode,x);
					if(fathnode == -1)
						continue;
					nofathernode = false;
					tailbridge += base[x];
					startnode = fathnode;
					break;
				}
			}
			for(int i = 0;i < numcolors;i++){
				msaseqs.push_back(tailbridge);
			}
//			cout << "tail bridge:" << tailbridge << endl;
			firstbub = false;
		}
		cout << "\nbubble id:" << bub +1<< endl;
		cout << "source/sink:" << source << ":" << sink << "\tpos:" << sourcepos << ":" << sinkpos << endl;
//step 1:Get the initial msa result
		cout << "step 1: Get the initial msa result" << endl;

//		int seqlength = sourcepos - sinkpos + kmer;
		string sourcelabel = dbg.node_label(source);
		char tmp = sourcelabel[kmer-1];
		string s(1,tmp);
		string emptystr(sinkpos - sourcepos,'-');
		string initallabel = s + emptystr;
		sourcelabel.append(emptystr);
		//string emptystr2(kmer+sinkpos-sourcepos,'-');
		string emptystr2(1+sinkpos-sourcepos,'-');
		vector<string> bubseq;
		vector<int> lengthcount(numcolors,0);
//		int maxlen = kmer + sinkpos - sourcepos;
		int maxlen = 1 + sinkpos - sourcepos;
//		cout << "startstr:" << sourcelabel << endl;
		for(int i = 0;i < numcolors;i++){
			if(bubcolor[i] == '1'){
				bubseq.push_back(initallabel);
				lengthcount[i]++;
			}
			else{
				bubseq.push_back(emptystr2);
			}
		}
//		vector<string> bubseq(numcolors,sourcelabel);
		vector<int> visited(dbg.num_nodes(),0);
		visited[source] = 1;
		queue<ssize_t> supernodelist;
		supernodelist.push(source);
		bubid.push_back(bub+1);
//		char* label_mat = new char[numcolors][maxlen];
		//char label_mat[maxlen * numcolors] = {'-'};
		char label_mat[maxlen * numcolors];
		for(int i = 0;i < maxlen * numcolors;i++){
			label_mat[i] = '-';
		}
//		char *label_mat = char label_mat[maxlen * numcolors]{'-'};
		//cerr << label_mat[2] << endl;
		for(int k = kmer-1;k < kmer;k++){
			char tmpchar = sourcelabel[k];
			for(int i = 0;i < numcolors;i++){
				if(researched_samples[i] == '0')
					continue;
				label_mat[i*maxlen] = tmpchar;
				//label_mat[i*maxlen+k] = tmpchar;
			}
		}
		while(!supernodelist.empty()){
			
			ssize_t startnode = supernodelist.front();
			ssize_t nextnode;
			for(unsigned long x = 1;x < dbg.sigma + 1;x++){
				if(startnode == sink){break;}
				nextnode = dbg.outgoing(startnode,x);
				if(nextnode == -1 || visited[nextnode] == 2)
					continue;
				string nodecolor = getnodecolor(nextnode);
				string commoncolor = ColorCap({nodecolor,bubcolor});
				if(commoncolor == zerocolor){continue;}
//			cout << startnode << "->" << nextnode << endl;
				if(visited[nextnode] == 0){
					//int nodelocate = nodepos[nextnode] -sourcepos + kmer - 1;
					int nodelocate = nodepos[nextnode] -sourcepos;
			//		char tmplabel = base[x];
					string nodelabel = dbg.node_label(nextnode);
					char tmplabel = nodelabel[kmer-1];	
			//		if(tmplabel != 'A' && tmplabel != 'C' && tmplabel != 'G' && tmplabel != 'T'){
			//			cout << "special char:" << tmplabel <<  endl;
			//		}
					int contain_sample = 0;
					for(int i = 0;i < numcolors;i++){
						if(commoncolor[i] == '1'){
							contain_sample++;
							string tmpseq = bubseq[i];
//							cout << nextnode << "\t" << nodepos[nextnode] << "\t" << i << endl;
							if(tmpseq[nodelocate] != '-' && tmpseq[nodelocate] != tmplabel){//meaning a sample has two nodes at the same position
								cout << "error: newbase is different from oldbase" << endl;
								cout << i << "\t" << nodelocate << endl;
								cout << tmpseq[nodelocate] << "\t" << tmplabel  << endl;
								cout << nextnode << "\tpos:" << nodepos[nextnode] << endl;
							}
							tmpseq[nodelocate] = tmplabel;
							bubseq[i] = tmpseq;
							label_mat[i*maxlen+nodelocate] = tmplabel;
							lengthcount[i]++;
						}
					}
					double rate = double(contain_sample)/double(numcolors);
					int tmppos = nodepos[nextnode];
					if(rate > 1){
						pos2keybase[tmppos] = tmplabel;
						pos2keybasevalue[tmppos] = rate;
			//		}else if(rate > 0.4){
			//			if((pos2keybase[tmppos] && rate > pos2keybasevalue[tmppos]) || !pos2keybase[tmppos]){
			//				pos2keybase[tmppos] = tmplabel;
			//				pos2keybasevalue[tmppos] = rate;
			//				
			//			}
					}
					supernodelist.push(nextnode);
					visited[nextnode] = 1;
				}
			} ///for(unsigned long x = 1;x < dbg.sigma + 1;x++)
			visited[startnode] = 2;
			supernodelist.pop();
			
		}//while(!supernodelist.empty()
//		int maxlength = Max(lengthcount);
		MyMSA << "bubble:" << bub+1 << endl;
		MyMSA << "source:" << source << "\tpos:" << sourcepos << "\tseq:" << dbg.node_label(source) << endl;;
		MyMSA << "sink:" << sink << "\tpos:" << sinkpos << "\tseq:" << dbg.node_label(sink) << endl;;
		for(int i = 0;i < numcolors;i++){
			MyMSA << bubseq[i] << endl;
		}
		bubblelens.push_back(lengthcount);
//		for(int i = 0;i < numcolors;i++){
//			for(int j = 0;j < maxlen;j++){
//				if(!label_mat[i*maxlen+j]){
//					cout << "-";
//				}else{
//					cout << label_mat[i*maxlen+j];
//				}
//			}
//			cout << endl;
//		}
//		cout << "Get initial msa results" << endl;


///////step 2:Start to prepare the data for fine-tune the result
		cout << "step 2: Start to prepare the data for fine-tune the result...." << endl;
	
		map<char,int> char2id = {{'A',1},{'C',2},{'G',3},{'T',4}};//- -> 0,A->1,C->2,G->3,T->4
		
		map<int,int> pos2vartype;
		//vector<int> gappos;
		//startpos = nodepos[source];
		bool existsgap = false;
//		cout << "initial position:" << endl;
		for(int i = 0;i < maxlen;i++){
			//int pos = sourcepos - kmer + 1 + i;
			int pos = sourcepos + i;
			vector<int> basecount(5,0);//五个位置分别对应 -,A,C,G,T
			for(int j = 0;j < numcolors;j++){
				if(researched_samples[j] == '0')
					continue;
				char tmpbase = label_mat[i + j * maxlen];
				if(!tmpbase||(tmpbase != 'A' && tmpbase!='C' && tmpbase!='G' && tmpbase!='T')){
					label_mat[i+j*maxlen] = '-';
					basecount[0]++;
				//	cout << "tmpbase: " << tmpbase << endl;
				}else{
					int tmpid = char2id[tmpbase];
					basecount[tmpid]++;
			//	cout << tmpbase << "\t" << tmpid << "\t" << basecount[tmpid] << endl;
				}
			}
			//cout << basecount[0] << "\t" << basecount[1] << "\t" << basecount[2] << endl;
			int snpcount = 0;
			for(int k = 1;k < 5;k++){
				if(basecount[k] > 0){
					snpcount++;
				}
			}
		//	cout << "basecount: " << basecount << endl;
			if(basecount[0] == 0){
				if(snpcount == 1){pos2vartype[pos] = 0;}
				else if(snpcount == 2){pos2vartype[pos] = 1;}
				else{pos2vartype[pos] = 2;}
//			cout << pos << "\t" << pos2vartype[pos] << endl;
			}
			else{
				if(snpcount == 1){pos2vartype[pos] = 3;}
				else{pos2vartype[pos] = 4;}
				existsgap = true;
//			cout << pos << "\t" << pos2vartype[pos] << endl;
		//		gapops.push_back(i);
			}
			if(snpcount == 0){cout << "get var type error at the position(raw): " << pos << endl;}
//			cout << pos << "\t" << pos2vartype[pos] << endl;
		}
//		for(int i = 0;i < numcolors;i++){
//			for(int j = 0;j < maxlen;j++){
//				cout << label_mat[i*maxlen+j];
//			}
//			cout << endl;
//		}
		
		//3.Determine whether we need to fine-tune or not
		//precondition: existsgap==true
		//input: pos2vartype //,gappos
		//output:wtregion (store the cutstart/cutend)
		cout << "step 3: Determine whether we need to fine-tune or not" << endl;
		cout << "if exists gap? " << existsgap << endl;
		if(existsgap == true){
			int gapstart = 0;
			int gapend = 0;
			bool ifwt = false;
			int gaplen = 0;
			int index = 0;
			vector<int> wtregion;
			while(index < maxlen){
				//int tmppos = sourcepos - kmer + 1 + index;
				int tmppos = sourcepos + index;
		//		cout << "index:" << index << endl;
				int vartype = pos2vartype[tmppos];
				//int vartype = pos2vartype[sourcepos + kmer - 1 - index];
				bool getstart = false;
				bool getend = false;
	//			cout << "index:" << index << "\tvar:" << vartype << endl;
				int type4 = 0;
				if(vartype == 3 || vartype == 4){
					if(gaplen == 0){
						gaplen++;
						gapstart = index;
					}
					else{
						gaplen++;
					}
					if(vartype == 4){type4++;}
					index++;
				}
				else{
					if(gaplen > 0){
			//			bool finetune = false;
						double rate = double(type4)/double(gaplen);
//						cout << "type4:"<< type4 << "\tgaplen:" << gaplen << "\trate:" << rate << endl;
						if(rate > 0.5){ifwt = true;}
						gapend = index;
						int twosidevarnum = 0;
						for(int i1 = 1;i1 <= 10;i1++){//这里10是人为设定
//							cout << "index: " << index << endl;
//							cout << tmppos+i1 << "\t" << pos2vartype[tmppos + i1] << endl;
					//		if(!pos2vartype[tmppos + i1]){
					//			cout << tmppos+i1 << " vartype doesnot exist" << endl;
					//			index = gapend + i1 -1;
					//			break;
					//		}
							if(pos2vartype[tmppos + i1] == 2 || pos2vartype[tmppos + i1] == 1){//snp情况大于2
								twosidevarnum++;
							}
							if(pos2vartype[tmppos + i1] > 2){
								twosidevarnum+=11;
								index = gapend + i1;
								break;
							}
						}
						if(twosidevarnum > 2 && twosidevarnum <= 10){
							ifwt = true;
						}
						else if(twosidevarnum <= 2 && gaplen < 3){
							index = gapend + 10;
						}
						else{
				//			gaplen++;
							ifwt = true;
						}
						gaplen = 0;
					//	if(gaplen > 1)
					//		ifwt = true;
						if(ifwt == true){
							int backlen = 0;
							int backpos = gapstart;
					//		cout << "backpos:" << backpos << endl;
							while(getstart == false){
								//if(pos2vartype[sourcepos -kmer+1 + backpos] > 1)
								if(pos2vartype[sourcepos + backpos] > 0){
									backlen = 0;
								}else{
									backlen++;
									if(backlen > 10){
										getstart = true;
									}
								
								}
								if(backpos < 1)
									break;
								backpos--;
							}
							gapstart = backpos;
					//		cout << "backpos:" << backpos << endl;
							int frontlen = 0;
							int frontpos = index;
					//		cout << "frontpos:" << frontpos << endl;
							while(getend == false){
								//if(pos2vartype[sourcepos -kmer+1 + frontpos] > 1)
								if(pos2vartype[sourcepos + frontpos] > 0){
									frontlen = 0;
								}else{
									frontlen++;
									if(frontlen > 10){
										getend = true;
									}
								
								}
								//if(frontpos == sinkpos - sourcepos + kmer - 1)
								if(frontpos == sinkpos - sourcepos)
									break;
								frontpos++;
							}
							gapend = frontpos-1;
					//		cout << "frontpos:" << frontpos << endl;
							if(gapend > maxlen-1){gapend = maxlen-1;}
							if(gapstart < 0){gapstart = 0;}
							int regionlen = wtregion.size();
							if(regionlen > 0 && gapstart < wtregion[regionlen - 1])
								gapstart = wtregion[regionlen - 1] + 1;
							wtregion.push_back(gapstart);
							wtregion.push_back(gapend);
							gapstart = 0;
							gapend = 0;
							index = frontpos;
				//			cout << "frontpos2:" << frontpos << endl;
						}
						
					}
					else{
						index++;
					}
				}
			}
			
			int wtsize = wtregion.size();
			if(wtsize > 0){
				cout << "need to fine-tune" << endl;
			}
			
			//4.prepare the data for fine-tune
			//input: wtregion, samp_labels
			//output:sortedseq(the first position is reference)
			//output:refhasgap(record that the reference has gap or not)
			cout << "step 4: prepare the data for fine-tune" << endl;

			ssize_t childnode = source;
			string longestpath = s;
			while(childnode != sink){
				ssize_t tmpstart = childnode;
				for(int k0 = 1;k0 < 5;k0++){
					ssize_t nextnode = dbg.outgoing(tmpstart,k0);
					if(nextnode == -1)
						continue;
					if(nodepos[nextnode] - nodepos[tmpstart] > 1)
						continue;
					string childlabel = dbg.node_label(nextnode);
					char tmpchar = childlabel[kmer-1];
					string s1(1,tmpchar);
					longestpath += s1;
					childnode = nextnode;
					break;
				}
			}
			cout << "longest path:" << longestpath << endl;

			vector<string> aligndsample_labels(numcolors,"");
			bool alignedresult = false;
			for(int i = wtsize/2-1;i >= 0;i--){
				int cutstart = wtregion[2 * i];
				int cutend = wtregion[2 * i + 1];
			
				cout << cutstart << "\t" << cutend << endl;
				map<string,vector<int>> rawseq2samples;
				map<string,vector<int>> seq2samples;
				map<string,int> seq2count;//the number of sequence
				map<string,int> seq2gaps;//the gap counts of the sequence
				vector<string> sortedseq;
				vector<string> alignedseq;//和sortedseq中的序列对应
				int nogap_num = 0;
				//int mincount = numcolors;
				//int maxcount = -1;
				vector<string> unsortseq;
				int seqcount = 0;//unique sequence counts
				int mingap = 999;
				for(int j = 0;j < numcolors;j++){//截取序列并做唯一化
					if(researched_samples[j] == '0')
						continue;
					string label = bubseq[j];
					string substr = label.substr(cutstart,cutend-cutstart+1);
//					cout << substr << endl;
					//string substr = label.substr(cutstart,cutstart-cutend+1);
					if(rawseq2samples.find(substr) != rawseq2samples.end()){
						vector<int> tmpsamp = rawseq2samples[substr];
						tmpsamp.push_back(j);
						rawseq2samples[substr] = tmpsamp;
						seq2count[substr]++;
					}
					else{
						vector<int> tmpsamp = {j};
						rawseq2samples[substr] = tmpsamp;
						seq2count[substr] = 1;
						int gapcount = 0;
						for(int k = 0;k <= cutend - cutstart + 1;k++){
							if(substr[k] == '-'){
					//		if(label_mat[j*maxlen+cutstart] == '-')
								gapcount++;
							}
						}
						seq2gaps[substr] = gapcount;
						if(gapcount == 0){nogap_num++;}
				//		cout << "mingap:" << mingap << ";" << gapcount << endl;
						if(mingap >= gapcount){ //to ensure the first element of unsortedseq has the minimum gap
							unsortseq.insert(unsortseq.begin(),substr);
							mingap = gapcount;
						}
						else{
							unsortseq.push_back(substr);
						}
						seqcount++;
					}
				}
				map<string,int> seq2sort;
				sortedseq.push_back(unsortseq[0]);
				seqcount = unsortseq.size();
		//	 	the priciple of sort: 1:less gap 2:more samples
				for(int j2 = 1;j2 < seqcount;j2++){
					int sortlen = sortedseq.size();
					string tmpstr = unsortseq[j2];
					int gapnum = seq2gaps[tmpstr];
					int seqcount = seq2count[tmpstr];
//					cout << j2 << "\t" << sortlen << "\t" << tmpstr << "\t" << gapnum << "\t" << seqcount << endl;
					bool addok = false;
	//					cout << endl;
	//					cout << gapnum << ":" << seqcount << "\t" << tmpstr << endl;
					for(int j3 = 0;j3 < sortlen;j3++){
						string tmpstr1 = sortedseq[j3];
						int gapnum1 = seq2gaps[tmpstr1];
						int seqcount1 = seq2count[tmpstr1];
		//				cout << gapnum1 << ":" << seqcount1 << "\t" << tmpstr1 << endl;
						if(gapnum1 > gapnum || (gapnum1 == gapnum && seqcount > seqcount1)){
							sortedseq.insert(sortedseq.begin()+j3,tmpstr);
							addok = true;
							break;
							
						}else{
							continue;
						}
					}
					if(addok == false){
						sortedseq.push_back(tmpstr);
					}
				}
//				string rawref = sortedseq[0];
		//		cout << "rawref:" << rawref << endl;
				for(int i0 = seqcount - 1;i0 >= 0;i0--){//去除序列中的‘-’
					string tmpstr = sortedseq[i0];
					vector<int> tmpsamples = rawseq2samples[tmpstr];
					int tmpcount = seq2count[tmpstr];
		//			cout << tmpcount << "\t" << tmpstr << endl;
					int tmpsize = tmpstr.size();
					for(int i1 = tmpsize-1;i1 >=0;i1--){
						if(tmpstr[i1] == '-'){
							tmpstr.erase(i1,1);
						}
					}
					sortedseq[i0] = tmpstr;
					if(seq2samples.find(tmpstr) != seq2samples.end()){
						sortedseq.erase(sortedseq.begin()+i0);
						vector<int> tmpsamp = seq2samples[tmpstr];
						int tmpsize = tmpsamp.size();
						for(int i1 = 0;i1 < tmpsize;i1++){
							int tmpsamp1 = tmpsamp[i1];
							if(find(tmpsamples.begin(),tmpsamples.end(),tmpsamp1) == tmpsamples.end()){
								tmpsamples.push_back(tmpsamp1);
							}
						}
					//	tmpsamp.insert(tmpsamp.end(),tmpsamples.begin(),tmpsamples.end());
					//	seq2samples[tmpstr] = tmpsamp;
					}
				//	else{
					seq2samples[tmpstr] = tmpsamples;
				//	}
					seq2count[tmpstr] = tmpcount;
				}
				//bool refhasgap = false;
//				cout << "nogap number:" << nogap_num << endl;
//
				map<int,char> tmppos2keybase;
				for(int i = cutstart;i <= cutend;i++){
					int tmppos = sourcepos + i;
					if(pos2keybase[tmppos] && !tmppos2keybase[i-cutstart-1]){
						cout << "keypos:" << tmppos << "\t" << pos2keybase[tmppos] << "\t" << i-cutstart << "\t" << cutstart << "\t" << endl;

						tmppos2keybase[i-cutstart] = pos2keybase[tmppos];
//						basescore[i-cutstart] = 1;
					}
				}
				//5.start to wt
				cout << "step 5: Start to fine-tune" << endl;
				vector<string> alignedseqs;
				string ref = longestpath.substr(cutstart,cutend-cutstart+1);
		//		string ref = sortedseq[0];
			//	cout << ref << endl;
				int reflen = ref.size();
				vector<int> basescore(reflen,1);
	//			map<int,char> tmppos2keybase;
			//	for(int i = cutstart;i <= cutend;i++){
			//		int tmppos = sourcepos + i;
			//		if(pos2keybase[tmppos]){
			//			cout << "keypos:" << tmppos << "\t" << pos2keybase[tmppos] << "\t" << i-cutstart << "\t" << cutstart << "\t" << reflen << endl;

//						tmppos2keybase[i-cutstart] = pos2keybase[tmppos];
			//			basescore[i-cutstart] = 1;
			//		}
			//	}
				seqcount = sortedseq.size();
				bool finishalign = false;
				int maxgap = 2;
				int multiple =  3;
			while(finishalign == false){
				for(int i = 0;i < seqcount;i++){
//					int seqlen = sortedseq[i].size();
					string seq = sortedseq[i];
					//	cout << "Before aligning1:\n" << ref << "\n" << sortedseq[i] << endl;
						Seqalign1_result aligned_result = SeqAlign1(seq,ref,tmppos2keybase,basescore,maxgap,multiple);
						vector<string> alignseqs = aligned_result.alignseq;
						tmppos2keybase = aligned_result.pos2keybase;
						basescore = aligned_result.basescore;
					//	cout << "Aligned result:" << endl;
					//	cout << alignseqs[0] << "\n" << alignseqs[1] << endl;
		//				cout << endl;
						string newref = alignseqs[0];
						seq = alignseqs[1];
						if(newref != ref){
							cout << "newref:\n" << newref << endl;
							ref = newref;
							alignedseqs.clear();
							break;
						}						
	//				cout << "basescore:" << basescore << endl;
					alignedseqs.push_back(seq);
				//	cout << seq << endl;
					if(i == seqcount-1)
						finishalign = true;
				}
			}
				//6.renew the result
				int seqnum = sortedseq.size();
		//		sortedseq[0] = rawref;
//				cout << "startpos:endpos\t" << sourcepos << ":" << sinkpos << "\t" << sourcepos-sinkpos+1 << endl;
//				cout << "cutstart:cutend\t" << cutstart << ":" << cutend << "\t" << cutend-cutstart+1 << endl;
				cout << "step 6: renew the result" << endl;
				int tmpcount = 0;
				for(int i0 = 0;i0 < seqnum;i0++){
					string rawseq = sortedseq[i0];
					
					string newseq = alignedseqs[i0];
				//	cout << rawseq << "\t" << newseq << endl;
					vector<int> samp = seq2samples[rawseq];
					int sampnum1 = samp.size();
					tmpcount += sampnum1;
					//cout << sampnum1 << "\t" << rawseq << "\t" << bubseq[samp[0]] << "\t" << newseq << endl;
//					cout << sampnum1 << "\t" << rawseq << "\t" << newseq << endl;
					for(int i1 = 0;i1 < sampnum1;i1++){
						int tmpsamp = samp[i1];
						string tmplabel = bubseq[tmpsamp];
						if(alignedresult != false){
				//			tmplabel = alignedsample_labels[tmpsamp];
						}
						tmplabel.erase(cutstart,cutend-cutstart+1);
						tmplabel.insert(tmplabel.begin()+cutstart,newseq.begin(),newseq.end());
						bubseq[tmpsamp] = tmplabel;
				//		cout << tmpsamp << ": " << tmplabel << endl;
					}
				} 
//				cout << "number:" << tmpcount << endl;
				//7.assess the result   alignedsample_labels
				cout << "step 7: assess the result" << endl;
				map<int,int> newpos2vartype;
				char newbase_mat[numcolors*(cutend-cutstart+1)];
				for(int i = 0;i < numcolors*(cutend-cutstart+1);i++){
					newbase_mat[i] = '-';
				}
				//int len = numcolors*(cutend-cutstart+1);
//				cout << len << endl;
				for(int i0 = 0;i0 < numcolors;i0++){
				//	string tmpstr = alignedsample_labels[i0];
					if(researched_samples[i0] == '0')
						continue;
					string tmpstr = bubseq[i0];
				//	cout << tmpstr << endl;
					for(int i1 = cutstart;i1 <= cutend;i1++){
						int tmppos = i0*(cutend-cutstart+1)+i1-cutstart;
				//		cout << tmppos << endl;
						newbase_mat[tmppos] = tmpstr[i1];
				//		cout << i0 << "\t" << i1 << "\t" << tmpstr[i1] << endl;
					}
				}
//				cout << "new base matrix:" << endl;
//				for(int i = 0;i < numcolors;i++){
//					for(int j = 0;j < cutend-cutstart+1;j++){
//						cout << newbase_mat[i*(cutend-cutstart+1)+j];
//					}
//					cout << endl;
//				}
//				cout << "cutstart->end " << cutstart << ":" << cutend << endl;
				vector<int> vartype_count(5,0);
		//		cout << "newpos:" << endl;
				for(int i3 = 0;i3 < cutend-cutstart+1;i3++){
					//int pos = sourcepos - kmer + 1 + cutstart + i3;
					int pos = sourcepos + cutstart + i3;
					vector<int> basecount(5,0);//五个位置分别对应 -,A,C,G,T
					for(int j = 0;j < numcolors;j++){
						if(researched_samples[j] == '0')
							continue;
						char tmpbase = newbase_mat[j*(cutend-cutstart+1)+i3];
						if(!tmpbase){
							newbase_mat[j*(cutend-cutstart+1)+i3] = '-';
							basecount[0]++;
						}else{
							int tmpid = char2id[tmpbase];
							basecount[tmpid]++;
						}
					}
					int snpcount = 0;
					for(int k = 1;k < 5;k++){
						if(basecount[k] > 0){
							snpcount++;
						}
					}
					if(basecount[0] == 0){//no gap,
						if(snpcount == 1){
							newpos2vartype[pos] = 0;
							vartype_count[0]++;
						}
						else if(snpcount == 2){
							newpos2vartype[pos] = 1;
							vartype_count[1]++;
						}
						else{
							newpos2vartype[pos] = 2;
							vartype_count[2]++;
						}
					}
					else{
						if(snpcount == 1){
							newpos2vartype[pos] = 3;
							vartype_count[3]++;
						}
						else{
							newpos2vartype[pos] = 4;
							vartype_count[4]++;
						}
						existsgap = true;
			//			gapops.push_back(i);
					}
					if(snpcount == 0){
						cout << "get var type error at the position(new)" << pos << endl;
						
					}//表示这一列全是gap
		//			cout << pos << "\t" << newpos2vartype[pos] << endl;
				}
				vector<int> rawvartype_count(5,0);
		//		cout << "rawpos:" << endl;;
				for(int i0 = cutstart;i0 <= cutend;i0++){
					//int pos = sourcepos - kmer + 1 + i0;
					int pos = sourcepos + i0;
					int tmpvartype = pos2vartype[pos];
		//			cout << pos << "\t" << tmpvartype << endl;;
					rawvartype_count[tmpvartype]++;
				}
				//然后比较rawvartype_count和vartype_count
				cout << "The improvement of fine-tune alignment:" << endl;
				cout << "vartype\t raw \t now" << endl;
				for(int i0 = 0;i0 < 5;i0++){//the larger,the more complex
					cout << i0 << "\t	" << rawvartype_count[i0] << "\t" << vartype_count[i0] << endl;
				}
			}
		}
		MyMSA << "fine-tune result" << endl;
		
		bubbleseqs.push_back(bubseq);
		string bridge = "";
		ssize_t startnode = source;
		bool nofathnode = false;
		int tmpx = -1;
		int branchedgenum = 0;
		while(dbg.indegree(startnode) == 1 && nofathnode == false){
			nofathnode = true;
//			cout << startnode << endl;
			string tmplabel = dbg.node_label(startnode);
			if(tmpx > -1){
				string s1(1,tmplabel[kmer-1]);
				bridge = s1 + bridge;
			}
			for(int x = 1;x < 5;x++){
				ssize_t fathnode = dbg.incoming(startnode,x);
				if(fathnode == -1)
					continue;
				startnode = fathnode;
				tmpx = x;
				nofathnode = false;
				branchedgenum++;
				break;
			}
			if(nofathnode == true){
				string head = tmplabel.substr(0,kmer-1);
				bridge = head + bridge;
			}
		}
//		cout << bubseq[numcolors - 1] << endl;
//		cout << bridge << endl;

		for(int i = 0;i < numcolors;i++){
		//	if(bridge == "" && branchedgenum == 0)
			if(pre_cutnode == true){
				int tmplen = bubseq[i].size();
				string newbubseq = bubseq[i].substr(0,tmplen-1);
				msaseqs[i] = newbubseq + msaseqs[i];
				msaseqs[i] = bridge + msaseqs[i];
			//	cout << newbubseq << endl;
			}else{
				msaseqs[i] = bubseq[i] + msaseqs[i];
				msaseqs[i] = bridge + msaseqs[i];
			}
			MyMSA << bubseq[i] << endl;
		}
		if(pre_cutnode == true)
			pre_cutnode = false;
		if(bridge == "" && branchedgenum == 0)
			pre_cutnode = true;
		//startbub = bub+1;
	}//for(int bub = 0;bub < setnum;bub++)
	mymsavar << "refpos\tnodepos\tref\talt\ttype\tvartype" << endl;
	int alignedlen = msaseqs[0].size();
	int refposcount = 0;
	int ts = 0;	
	int sp = 0;
	for(int i = 0;i < alignedlen;i++){
		char refbase = '*';
		vector<int> basecount(5,0);
		for(int j = 0;j < numcolors;j++){
			string sampstr = msaseqs[j];
			char tmpbase = sampstr[i];
			if(tmpbase == 'A'){
				basecount[0]++;
			}else if(tmpbase == 'C'){
				basecount[1]++;
			}else if(tmpbase == 'G'){
				basecount[2]++;
			}else if(tmpbase == 'T'){
				basecount[3]++;
			}else{
				basecount[4]++;
			}
			if(j == numcolors - 1){
				refbase = tmpbase;
				if(tmpbase != '-')
					refposcount++;
			}
		}
		//cout << i+1 << "\t" << basecount << endl;
		if(basecount[0] != numcolors && basecount[1] != numcolors && basecount[2] != numcolors && basecount[3] != numcolors){
		//	cout << i+1 << "\t" << basecount << endl;
			if(refbase == '-'){
				mymsavar << refposcount+1 << "\t" << i+1  << "\t" << refbase << "\t";
			}else{
				mymsavar << refposcount << "\t" << i+1  << "\t" << refbase << "\t";
			}
			int tmpvartype = -1;//1-3 -> 2-4 type bases
			int tmpgapid = 4;
			for(int j1 = 0;j1 < 4;j1++){
				if(basecount[j1]== 0){
					continue;
				}
				if(base[j1+1] == refbase){
					tmpvartype++;
					continue;
				}
				mymsavar << base[j1+1];
				tmpgapid = j1+1;
				tmpvartype++;
				break;
			}
			for(int j1 = tmpgapid;j1 < 4;j1++){
				if(basecount[j1]== 0){
					continue;
				}
				if(base[j1+1] == refbase){
					tmpvartype++;
					continue;
				}
				mymsavar << ":" << base[j1+1];
				tmpvartype++;
			}
			if(refbase == '-'){
				tmpvartype += 8;
			}else{
				if(basecount[4] > 0){
					tmpvartype+=4;
				}
			}
			mymsavar << "\t" << tmpvartype << "\t";
			if(tmpvartype < 4){
				mymsavar << "SNP" << endl;
			}else if(tmpvartype < 8){
				mymsavar << "Del" << endl;
			}else{
			 	mymsavar << "Ins" << endl;
			}
		}else{
			ts++;
		}
		/// caculate TS score and SP score
		/// basecount
		for(int i1 = 0;i1 < 4;i1++){
			if(basecount[i1] == 0)
				continue;
			sp = sp + (basecount[i1]*(basecount[i1]-1))/2;
//			for(int i2 = i1+1;i2 < 5;i2++){
//				if(basecount[i2] == 0)
//					continue;
//				sp = sp - basecount[i1]*basecount[i2];
//			}
		}
		
	}
	double TC = double(ts)/double(alignedlen);
	double SP = double(sp)/double(double(alignedlen)*double(numcolors)*(double(numcolors)-1)/2);
	cout << "Same colunm:" << ts << endl;
	cout << "Aligned length:" << alignedlen << endl;
	cout << "Total Column(TC):" << TC << endl;
	cout << "Same pairs:" << sp << endl;
	cout << "Sum of Pairs(SP):" << SP << endl;
	for(int i = 0;i < numcolors;i++){
		mymergeMSA << msaseqs[i] << endl;
	}
	cout << "MSA are listed in the myMSA.txt\n";
	cout << refposcount << endl;
	msa_set.bubid = bubid;
	msa_set.bubbleseqs = bubbleseqs;
	msa_set.bubblelens = bubblelens;
	int totaltime = getMilliSpan(t);
	if(totaltime < 1000){
		cout << "Get MSA's time: " << getMilliSpan(t) << "ms" << endl;
	}
	else{
		int tolsecond = totaltime/1000;
		int hour = tolsecond/3600;
		int second1 = tolsecond%3600;
		int minute = second1/60;
		int second = second1%60;
		cout << "Get MSA's time: " << hour << "h " << minute << "m " << second << "s " << endl;
	}
	return msa_set;

}


//####################################################################################################################################################
//###########################################################  Variants Analysis  ####################################################################
//####################################################################################################################################################

void variants_analysis(MSA msa){
	vector<int> bubids = msa.bubid;
	vector<vector<string>> bubbleseqs = msa.bubbleseqs;
	vector<vector<int>> bubblelens = msa.bubblelens;
	vector<vector<ssize_t>> pairnodesets = msa.pairnodesets;
	map<ssize_t,int> nodepos = msa.nodepos;
	vector<string> bubblecolors = msa.bubblecolors;
	int kmer = dbg.k - 1;
	int numcolors = colors.size()/dbg.size();
	string stdcolor(numcolors,'1');
	int bubblenum = bubids.size();
	vector<vector<string>> alter_vars;
	vector<int> vars_pos;
	vector<int> typevarcount(2,0);
	int maxvarlen = -1;
//	int type1var = 0;
//	int type2var = 0;
	cout << "Start to find variants..." << endl;
	cout << "newid\toldid\tposition\trefvar\tvariation..." << endl;
	int totalvar = 0;
	int totalvarbase = 0;
	vector<int> totalvartypecount(3,0);
//	int skip = 3;
	for(int bub = 0;bub < bubblenum;bub++){
		vector<string> bubbleseq = bubbleseqs[bub+1];
		vector<int> bubblelen = bubblelens[bub+1];
//		int maxlen = Max(bubblelen);
	//	string longestsamp = stdcolor;
	//	if(lengthcount[i] < maxlen){
	//		longestsamp[i] = 0;
	//	}
		int bubid = bubids[bub+1];
		vector<ssize_t> pairnodeset = pairnodesets[bubid];
		ssize_t source = pairnodeset[1];
	//	ssize_t sink = pairnodeset[2];
		int sourcepos = nodepos[source];
//		int sinkpos = nodepos[sink];
		int row = bubbleseq.size();
		int col = bubbleseq[0].size();
//		cout << bubid << "\t" << bubblelen << "\t" << col << "\tsource:" << source << endl;
		int ord = 0;
		while(ord < col){
			int i = 0;//record the visiting sample's id
			int skip = 1;//the length of variant at this position
			map<string,int> varcount;
			bool meetgap = false;
			string refvar;
			while(i < row){
				if(researched_samples[i] == '0'){
					i++;
					continue;
				}
				string bubseq = bubbleseq[i];
				if(meetgap == false && bubseq[ord] == '-'){
					i = 0;
					varcount.clear();
					meetgap = true;
					continue;
				}
				if(meetgap == true && bubseq[ord+skip] == '-'){//there are more contiguous gaps
					skip++;
					i = 0;
					varcount.clear();
					continue;
				}
				string varbase = bubseq.substr(ord,skip);
				if(varcount[varbase] == false){
					varcount[varbase] = 1;
				}
				else{
					varcount[varbase]++;
				}
				i++;
				if(i == row){
					refvar = varbase;
				}
			}
			int varnum = varcount.size();
			int varpos = sourcepos + ord - kmer + 1;
			ord += skip;
			if(skip > maxvarlen){
				maxvarlen = skip;
			}
//			cout << "There are " << varnum << " variants" << endl;
			if(varnum >= 2){
//				cout << bubid << " has more than two variants at the position " << varpos << endl;
				totalvar++;
				int typevarcountsize = typevarcount.size();
				if(varnum-1 > typevarcountsize){
					for(int i0 = 0;i0 < varnum - 1 - typevarcountsize;i0++){
						typevarcount.push_back(0);
					}
				}
				typevarcount[varnum-2]++;
				map<string,int>::iterator iter;
				iter = varcount.begin();
	//			auto iter = varcount.begin();
				vector<string> alter_var;
				int varlen = refvar.size();
				totalvarbase+=varlen;
				vector<int> vartype(varlen,0);
				for(int i = 0;i < varlen;i++){
					if(refvar[i] == '-')
						vartype[i] = 2;
				}
				cout << bub+1 << "\t" << bubid << "\t" << varpos << "\t" << refvar << "(" << varcount[refvar] << ")";
				while(iter != varcount.end()){
			//	for(auto iter = varcount.begin();iter < varcount.end();iter++){
					string alter = iter -> first;
					int altercount = iter -> second;
					alter_var.push_back(alter);
					cout << "\t" << alter << "(" << altercount << ")";
					iter++;
					for(int i = 0;i < varlen;i++){
						if(alter[i] != refvar[i] && vartype[i] != 2 && alter[i] == '-')
							vartype[i] = 2;
						if(alter[i] != refvar[i] && vartype[i] == 0 && alter[i] != '-')
							vartype[i] = 1;
					}
				}
				vector<int> vartypecount(3,0);
				for(int i = 0;i < varlen;i++){
					int tmptype = vartype[i];
					vartypecount[tmptype]++;
					totalvartypecount[tmptype]++;
				}
				alter_vars.push_back(alter_var);
				cout << endl;
			}
			else if(varnum == 1){
				//means that there is not variants at this position
			}
			else{cout << bubid << " has error at the position " << varpos << endl;}
		}
	}
	cout << "There are " << totalvar << " variants" << endl;
	int maxvarcount = typevarcount.size();
	cout << "genetype counts(2-" << maxvarcount + 1 << "):" << typevarcount << endl;
	cout << "max variation length: " << maxvarlen << endl;
	cout << "Total variation base count: " << totalvarbase << endl;
	cout << "Total variation type count: SNP(" << totalvartypecount[1] << ")\tIndel(" << totalvartypecount[2] << ")" << endl;
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

	int numcolors = colors.size()/dbg.size();
	string zerocolor(numcolors,'0');
	researched_samples.assign(numcolors,'1');
//	string researched_samples_file = "researched_samples.txt";
//	researched_samples = read_researched_samples_from_file(researched_samples_file);
	cout << researched_samples << endl;
//	nearby();

	C_Su c_su = read_mypairs_from_file("cSupB_topo.txt");
//	vector<vector<ssize_t>> pairnodeset = c_su.pairnodesets;
//	vector<string> bubblecolor = c_su.bubblecolors;
//	vector<ssize_t> pairnodes = pairnodeset[1];
//	cout << "pairnodes:" << pairnodes[0] << "\t" << pairnodes[1] << ":" << pairnodes[2] << "\t" << bubblecolor[1] << endl;
	
	MSA msa = GetMSA(c_su);
//	variants_analysis(msa);
//	Cytoscape();

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
