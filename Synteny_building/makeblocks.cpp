#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <regex>
#include <dirent.h>
#include <unistd.h>
#include <list>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <limits.h>

#define MINOVL 0.4
#define AFEW 0.3
#define MINDESSEG 0.05 // min 
#define MINOUTSEG 0.02 // min
#define MAXDEP 30 //max level for reading net file 

using namespace std;

vector < string > spc_list;
vector < int > spctag_list; // 
vector < string > chr_name;

int min_rsl = -1, ref_idx = 0, spc_cnt = 0;
int maxchr = 0;
int maxorder = 0;
string chDir, netDir, ref_spc;
string outDir;
vector < string > raw_segs_onmem;

enum segstate {
    FIRST = 0, LAST, BOTH, MIDDLE
};
struct seg {
    char *aChr, *bChr;
    int aBeg, aEnd, bBeg, bEnd;
    int chainID;
    char dir;
    struct seg *next;
};
void free_seg(seg* head_ref)  {        
    seg* current = head_ref;  
    seg* next;  
  
    while (current != NULL)  {  
        next = current->next;  
        free(current->aChr);
        free(current->bChr);
        delete(current);  
        current = next;  
    } 
}
struct block {
    char *refchrom;
    int refbeg, refend;
    struct seg *(*speseg);
    struct block *next;
};
void free_block(block* head_ref)  {        
    block* current = head_ref;  
    block* next;  
  
    while (current != NULL)  {  
        next = current->next;
        for(int i = 0 ; i < spc_cnt ; i ++){
            free_seg(current->speseg[i]);
        }
        delete(current);  
        current = next;  
    } 
}

struct gf_list {
	int size, fgap, sgap;
	struct gf_list *next;
};

struct chain_list {
	int cid;
	char *fchr, *schr;
    int fbeg, fend, sbeg, send, slen;
	char forient, sorient;
	struct gf_list *gf;
	struct chain_list *next;
};

map < int , struct chain_list* > chain_list;
map < int, string > hs_exist;

struct block *my_allocate_newblock(){
    struct block *newblock;
    int i ;
    newblock = new block;
    newblock->next = NULL;
    newblock->speseg = (struct seg **)malloc(sizeof(struct seg*) * spc_cnt);
    for (i = 0 ; i < spc_cnt ; i ++){
        newblock->speseg[i] = NULL;
    }
    newblock->refbeg = INT_MAX;
    newblock->refend = INT_MIN;
    newblock->refchrom = nullptr;
    return newblock;
}
struct seg_list {
    int id, beg, end, subid, chid, chnum;
    int *cidlist;
    char *chr;
    char orient;
    enum segstate state;
    struct seg_list *next;
};
void free_seg_list (seg_list* head_ref)  {        
    seg_list* current = head_ref;  
    seg_list* next;  
  
    while (current != NULL)  {  
        next = current->next;  
        if (current->cidlist != NULL){
            free(current->cidlist);
        }
        free(current->chr);  
        delete(current);  
        current = next;  
    } 
}
struct block_list {
    int id, isdup;
    int left, right;
    struct seg_list *(*speseg);
    struct block_list *next;
};
struct block_list *my_allocate_newblock_list(){
    
    struct block_list *newblock;
    int i ;
    newblock = new block_list;
    newblock->next = NULL;
    newblock->speseg = (struct seg_list **)malloc(sizeof(struct seg_list*) * spc_cnt);
    for (i = 0 ; i < spc_cnt ; i ++){
        newblock->speseg[i] = NULL;
    }
    newblock->isdup = 0;
    newblock->id=0;
    return newblock;
}
void free_block_list(struct block_list *blk) {
    struct block_list *p, *q;
    int i;
    p = blk;
    for (;;) {
        q = p->next;
        for (i = 0; i < spc_cnt; i++)
            if (p->speseg[i] != NULL)
                free_seg_list(p->speseg[i]);
        delete(p);
        if (q == NULL)
            break;
        else
            p = q;
    }
}


string CHOMP(string str) {
    string::size_type pos = str.find_last_not_of("\n \t");
    if (pos != string::npos)
        str = str.substr(0, pos + 1);
    return str;
}
int STR2INT(string s) {
    return atoi(s.c_str());
}
int get_spc_idx(string spcname) {
    int i = 0;
    vector < string > ::iterator it = find(spc_list.begin(), spc_list.end(), spcname);
    i = distance(spc_list.begin(), it);
    return i;
}
vector < string > SPLIT(string s, char delimiter = '1') {
    vector < string > return_arr;
    istringstream iss(s);
    if (delimiter == '1') { //default, split with spaces
        do {
            string sub;
            iss >> sub;
            return_arr.push_back(sub);
        } while (iss);

        return_arr.pop_back();
        return return_arr;
    } else {
        string temp;
        while (getline(iss, temp, delimiter)) {
            return_arr.push_back(temp);
        }
        return return_arr;
    }
}
void READ_CONFIG(string f_config) {
    int tag;
    char sn[1000];
    char buf[PATH_MAX];
    FILE * configFile = fopen(f_config.c_str(), "r");
    int idx = 0;
    int spc_index = 0;
    while (fgets(buf, PATH_MAX, configFile)) {
        if (buf[0] == '\n' || buf[0] == '#') {
            continue;
        }
        if (buf[0] == '>' && strstr(buf, "netdir") != NULL) {
            fgets(buf, PATH_MAX, configFile);
            netDir = CHOMP(string(buf));
        } else if (buf[0] == '>' && strstr(buf, "chaindir") != NULL) {
            fgets(buf, PATH_MAX, configFile);
            chDir = CHOMP(string(buf));
        } else if (buf[0] == '>' && strstr(buf, "resolution") != NULL) {
            fgets(buf, PATH_MAX, configFile);
            min_rsl = STR2INT(string(buf));
        } else if (buf[0] == '>' && strstr(buf, "species") != NULL) {
            while (fgets(buf, PATH_MAX, configFile)) {
                if (buf[0] == '\n') {
                    break;
                }
                if (sscanf(buf, "%s %d", sn, & tag) != 2) {
                    break;
                }
                spc_list.push_back(string(sn));
                spctag_list.push_back(tag);
                if (tag == 0) {
                    ref_spc = string(sn);
                    ref_idx = idx;
                }
                idx++;
            }
        }
    }
    spc_cnt = spc_list.size();
    if (spc_list.size() != spctag_list.size()) {
        fprintf(stderr, "[ERROR] Species list is not properly prepared\n");
        exit(0);
    }
} //read configs
int GET_LEVEL(string s) {
    int i = 0;
    while (s.at(i) == ' ') {
        i++;
        if (i > 30) {
            fprintf(stderr, "[ERROR] MAXDEP = %d\n", 30);
        }
    }
    return i;
}
void READ_NETS(int spcidx, string netdir) {
    DIR * dir;
    struct dirent * ent;
    int fbeg = 0, flen = 0, sbeg = 0, slen = 0, cid = 0, level = 0;
    string chrom, ref_chr;
    char orient;
    map < int, int > val;
    map < int, int > fgapbeg, fgapend, sgapbeg, sgapend;
    map < int, string > gapchrom;
    map < int, char > gaporient;
    vector < string > tmp_arr;
    string this_raw_seg = "";

    string raw_segs = outDir + "/" + spc_list[spcidx] + ".raw.segs";
    FILE * fp1 = fopen(raw_segs.c_str(), "w");
    vector <string> filelist ;
    if ((dir = opendir(netdir.c_str())) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            string netName = ent -> d_name;
            if (netName == "." || netName == "..") {
                continue;
            }
			string formatCheck = netName.substr(netName.length()-4,netName.length());
			if (formatCheck != ".net"){
				continue;
			}
            string ref_chr = netName.substr(0, netName.length() - 4); // get rid of ".net"
            string netFile = netdir + "/" + netName;
            filelist.push_back(netFile);
        } //while
        closedir(dir);
    } else {
        fprintf(stderr, "[Error] Fail to open net directory \n");
        exit(0);
    }
    sort(filelist.begin(), filelist.end());
    for(int kk = 0 ; kk < filelist.size(); kk ++){
            string netFile = filelist[kk];
            char buf[5000];
            FILE * netstream = fopen(netFile.c_str(), "r");
            for (int i = 0; i < MAXDEP; i++) {
                val[i] = 0;
            }
           
            while (fgets(buf, 5000, netstream)) {
                if (buf[0] == '\n' || buf[0] == '#') {
                    continue;
                }
                if (buf[0] == 'n') { // net 
                    char refchrom[500];
                    if (sscanf(buf, "net %s %*s", refchrom) != 1) {
                        fprintf(stderr, "cannot parse: %s", buf);
                        exit(1);
                    }
                    ref_chr = string(refchrom);
                    chr_name.push_back(ref_chr);
                }
                char type[10];
                if (sscanf(buf, "%s %*s", type) != 1) {
                    fprintf(stderr, "cannot parse: %s", buf);
                    exit(1);
                }
                if (type[0] == 'f') { // fill 
                    level = GET_LEVEL(string(buf));
                    level /= 2;
                    for (int i = level; i < MAXDEP; i++) {
                        val[i] = 0;
                    }
                    char fchrom[500];
                    if (sscanf(buf, "%*s %d %d %s %c %d %d id %d %*s", & fbeg, & flen, fchrom, & orient, & sbeg, & slen, & cid) != 7) {
                        fprintf(stderr, "cannot parse: %s", buf);
                        exit(1);
                    }
                    chrom = string(fchrom);
                    if (flen > min_rsl || slen > min_rsl) {
                        val[level] = 1;

                        fprintf(fp1, "%d s %s.%s:%d-%d %s.%s:%d-%d %c %d", level, ref_spc.c_str(), ref_chr.c_str(), fbeg, fbeg + flen, spc_list[spcidx].c_str(), chrom.c_str(), sbeg, sbeg + slen, orient, cid);
                        this_raw_seg = "";
                        this_raw_seg = to_string(level) + " s " + ref_spc + ":" + ref_chr + ":" + to_string(fbeg) + "-" + to_string(fbeg + flen) + " " + spc_list[spcidx] + ":" + chrom + ":" + to_string(sbeg) + "-" + to_string(sbeg + slen) + " " + orient + " " + to_string(cid);
                        if (level == 0) {
                            raw_segs_onmem.push_back(this_raw_seg);
                            fprintf(fp1, "\n");
                        } else {
                            int k;
                            for (k = level - 1; k >= 0; k--) {
                                if (val[k] == 1) {
                                    break;
                                }
                            }
                            if (k < 0) {
                                fprintf(fp1, " [NP]\n");
                                this_raw_seg += " [NP]";
                                raw_segs_onmem.push_back(this_raw_seg);
                            } else {
                                fprintf(fp1, " [%d %d %s %d %d %c]\n", fgapbeg[k], fgapend[k], gapchrom[k].c_str(), sgapbeg[k], sgapend[k], gaporient[k]);
                                this_raw_seg += " [ " + to_string(fgapbeg[k]) + " " + to_string(fgapend[k]) + " " + gapchrom[k] + " " + to_string(sgapbeg[k]) + " " + to_string(sgapend[k]) + " " + gaporient[k] + " ]";
                                raw_segs_onmem.push_back(this_raw_seg);
                            }
                        }
                    }
                } else if (type[0] == 'g') {
                    level = GET_LEVEL(string(buf));
                    level /= 2;
                    --level;
                    char ar_gapchrom[500];
                    if (sscanf(buf, "%*s %d %d %s %c %d %d %*s", & (fgapbeg[level]), & (fgapend[level]), ar_gapchrom, & (gaporient[level]), & (sgapbeg[level]), & (sgapend[level])) != 6) {
                        fprintf(stderr, "cannot parse: %s", buf);
                        exit(1);
                    }
                    gapchrom[level] = string(ar_gapchrom);
                    fgapend[level] += fgapbeg[level];
                    sgapend[level] += sgapbeg[level];
                    if (sgapend[level] - sgapbeg[level] > min_rsl) {
                        this_raw_seg = "";
                        this_raw_seg = to_string(level) + " g " + ref_spc + ":" + ref_chr + ":" + to_string(fgapbeg[level]) + "-" + to_string(fgapend[level]) + " " + spc_list[spcidx] + ":" + gapchrom[level] + ":" + to_string(sgapbeg[level]) + "-" + to_string(sgapend[level]) + " " + gaporient[level];
                        fprintf(fp1, "%d g %s.%s:%d-%d %s.%s:%d-%d %c\n", level, ref_spc.c_str(), ref_chr.c_str(), fgapbeg[level], fgapend[level], spc_list[spcidx].c_str(), gapchrom[level].c_str(), sgapbeg[level], sgapend[level], gaporient[level]);
                        raw_segs_onmem.push_back(this_raw_seg);
                    }
                }

            } //while for a file 
            fclose(netstream);        
    }
    fclose(fp1);
}
struct seg * PROC_SEGS(int spcidx){
    struct seg *slist, *p, *q, *r, *prp , *plast;
    slist = plast = p = r = prp = NULL;
    int tobreak = 0, level = 0;
    string type, fchrom, schrom;
    int fgbeg, fgend, sgbeg, sgend;
    int b1, e1, b2, e2;
    string gchrom;
    char gorient;

    vector < string > tmp_arr;
    regex reg_raw("(\\w+):(\\S+):(\\d+)-(\\d+)");
    smatch reg_arr;

    for (int i = 0; i < raw_segs_onmem.size(); i++) {
        tmp_arr = SPLIT(raw_segs_onmem[i]);
        level = STR2INT(tmp_arr[0]);
        type = tmp_arr[1];
        if (level == 0 && type == "s") {
            p = new(struct seg);
            p->dir = (tmp_arr[4])[0];
            p->next = NULL;
            p->chainID = STR2INT(tmp_arr[5]);
            regex_search(tmp_arr[2], reg_arr, reg_raw);
            p->aChr = (char*)malloc(sizeof(char)*reg_arr[2].length()+1);
            strcpy(p->aChr,reg_arr[2].str().c_str());
            p->aBeg = STR2INT(reg_arr[3]);
            p->aEnd = STR2INT(reg_arr[4]);
            regex_search(tmp_arr[3], reg_arr, reg_raw);
            p->bChr = (char*)malloc(sizeof(char)*reg_arr[2].length()+1);
            strcpy(p->bChr,reg_arr[2].str().c_str());
            p->bBeg = STR2INT(reg_arr[3]);
            p->bEnd = STR2INT(reg_arr[4]);
            if (slist == NULL){
                slist = plast = p;
            }else{
                while(plast->next != NULL){
                    plast = plast->next;
                }
                plast->next = p;
                plast = p;
            }
        }else{
            if(type == "s"){
                p = new(struct seg);
                p->dir = (tmp_arr[4])[0];
                p->chainID = STR2INT(tmp_arr[5]);
                regex_search(tmp_arr[2], reg_arr, reg_raw);
                p->aChr = (char*)malloc(sizeof(char)*reg_arr[2].length()+1);
                strcpy(p->aChr,reg_arr[2].str().c_str());
                p->aBeg = STR2INT(reg_arr[3]);
                p->aEnd = STR2INT(reg_arr[4]);
                regex_search(tmp_arr[3], reg_arr, reg_raw);
                p->bChr = (char*)malloc(sizeof(char)*reg_arr[2].length()+1);
                strcpy(p->bChr,reg_arr[2].str().c_str());
                p->bBeg = STR2INT(reg_arr[3]);
                p->bEnd = STR2INT(reg_arr[4]);
                if(tmp_arr[6]  != "[NP]") {
                    fgbeg = STR2INT(tmp_arr[7]);
                    fgend = STR2INT(tmp_arr[8]);
                    gchrom = tmp_arr[9];
                    sgbeg = STR2INT(tmp_arr[10]);
                    sgend = STR2INT(tmp_arr[11]);
                    gorient = tmp_arr[12].c_str()[0];
                }
                tobreak = 0;
                for (q = plast, prp = NULL ; q != NULL; q = q->next) {
                    if(strcmp(p->aChr,q->aChr)==0 && p->aBeg >= q->aBeg && p->aEnd <= q->aEnd){
                        tobreak = 1;
                        break;
                    }
                    if(prp != NULL && strcmp(p->aChr,q->aChr)==0 && p->aBeg >= prp->aEnd && p->aEnd <= q->aBeg){
                        tobreak=0;
                        break;
                    }
                    prp = q;
                }
                if(tobreak ==1){
                    r = new(struct seg);
                    p->next =r;
                    r->next =q->next;
                    q->next =p;

                    r->aEnd = q->aEnd;
                    q->aEnd = fgbeg;
                    r->aBeg = fgend;
                    r->aChr = (char*)malloc(sizeof(char)*string(q->aChr).length()+1);
                    strcpy(r->aChr,q->aChr);
                    r->bChr = (char*)malloc(sizeof(char)*string(q->bChr).length()+1);
                    strcpy(r->bChr,q->bChr);
                    r->chainID = q->chainID;
                    r->dir = q->dir;
                    if (q->dir == '+') {
                        r->bEnd = q->bEnd;
                        q->bEnd = sgbeg;
                        r->bBeg = sgend;
                    } else {
                        r->bBeg = q->bBeg;
                        q->bBeg = sgend;
                        r->bEnd = sgbeg;
                    }
                }else{ // tobreak == 0;
                    if(prp){
                        p->next = prp->next;
                        prp->next = p;
                    }else{
                        prp = p;
                    }
                }
            }else{//type = "g"
                tobreak = 0;
                regex_search(tmp_arr[2], reg_arr, reg_raw);
                fchrom = reg_arr[2];
                b1 = STR2INT(reg_arr[3]);
                e1 = STR2INT(reg_arr[4]);
                regex_search(tmp_arr[3], reg_arr, reg_raw);
                schrom = reg_arr[2];
                b2 = STR2INT(reg_arr[3]);
                e2 = STR2INT(reg_arr[4]);
                for(q = plast, prp = NULL; q != NULL; q = q->next){
                    if (string(q->aChr) == fchrom && b1 >= q->aBeg && e1 <= q->aEnd){
                        tobreak=1;
                        break;
                    }
                    if(prp != NULL && fchrom == string(q->aChr) && e1 <= q->aBeg && b1 >= prp->aEnd){
                        tobreak=0;
                        break;
                    }         
                    prp = q; 
                }
                
                if(tobreak == 1 && q->bChr == schrom){
                    r = new(struct seg);
                    r->next = q->next;
                    q->next = r;
                    r->aEnd = q->aEnd;
                    q->aEnd = b1;
                    r->aBeg = e1;
                    r->aChr = (char*)malloc(sizeof(char)*string(q->aChr).length()+1);
                    strcpy(r->aChr,q->aChr);
                    r->bChr = (char*)malloc(sizeof(char)*string(q->bChr).length()+1);
                    strcpy(r->bChr,q->bChr);
                    r->chainID = q->chainID;
                    r->dir = q->dir;
                    if (q->dir == '+') {
                        r->bEnd = q->bEnd;
                        q->bEnd = b2;
                        r->bBeg = e2;
                    } else {
                        r->bBeg = q->bBeg;
                        q->bBeg = e2;
                        r->bEnd = b2;
                    }
                }
            }
        }
    }
    string proc_segs = outDir + "/" + spc_list[spcidx] + ".processed.segs";
    FILE * fp2 = fopen(proc_segs.c_str(), "w");
    string prv_chr = "";
    for(p = slist; p !=NULL ; p = p->next){
        if((string(p->aChr) != prv_chr)){
            fprintf(fp2,"#\n");
        }
        prv_chr = string(p->aChr);
        fprintf(fp2, "%s.%s:%d-%d %s.%s:%d-%d %c %d\n", ref_spc.c_str(), p->aChr, p->aBeg, p->aEnd, spc_list[spcidx].c_str(), p->bChr, p->bBeg, p->bEnd, p->dir, p->chainID);
    }
    fclose(fp2);
    raw_segs_onmem.clear();
    return slist;
}
void fill_block(struct block *blck, int idx, struct seg *sg){
    if (!(blck->refchrom)){
        blck->refchrom = (char*)malloc(sizeof(char)*string(sg->aChr).length()+1);
        strcpy(blck->refchrom,sg->aChr);
    }else if (strcmp(blck->refchrom,sg->aChr)!=0){
        fprintf(stderr,"CHROM DISAGREE: %s %s", blck->refchrom, sg->aChr);
    }
    
    blck->refbeg = min(blck->refbeg, sg->aBeg);
    blck->refend = max(blck->refend, sg->aEnd);
    struct seg *newsg = new seg;
    newsg->next = NULL;
    newsg->aChr = (char*)malloc(sizeof(char)*string(sg->aChr).length()+1);
    newsg->bChr = (char*)malloc(sizeof(char)*string(sg->bChr).length()+1);
    strcpy(newsg->aChr, sg->aChr);
    strcpy(newsg->bChr, sg->bChr);
    newsg->aBeg = sg->aBeg;
    newsg->aEnd = sg->aEnd;
    newsg->bBeg = sg->bBeg;
    newsg->bEnd = sg->bEnd;
    newsg->dir = sg->dir;
    newsg->chainID = sg->chainID;

    blck->speseg[idx] = newsg;
}
void find_insert_position(struct seg *sg, struct block *blockhead, struct block *last, struct block **prv, struct block **nxt, struct block **fst, struct block **lst){
    struct block *p, *pp;
    *prv = *nxt = *fst = *lst = NULL;
    pp = NULL;
    if (last != NULL)
        p = last;
    else
        p = blockhead;
    for (; p != NULL; p = p->next) {
        if (strcmp(p->refchrom,sg->aChr)!=0 && p->next != NULL && strcmp(p->next->refchrom,sg->aChr)==0){
            pp = p;
        }else if(strcmp(p->refchrom,sg->aChr)==0){
            if((p->next != NULL && strcmp(p->next->refchrom,sg->aChr)==0 && p->refend <= sg->aBeg && sg->aEnd <= p->next->refbeg) || ((p->next == NULL || strcmp(p->next->refchrom,sg->aChr)!=0) &&  p->refend <= sg->aBeg)){
                *prv = p;
                *nxt = p->next;
                *fst = *lst = NULL;   
                break;
            }
            if(((pp != NULL && p == pp->next) || (pp == NULL && p == blockhead)) &&  p->refbeg >= sg->aEnd){
                *prv = pp;
                *nxt = p;
                *fst = *lst = NULL; 
                break;
            }
            // determine prv and fst
            if (p->refbeg <= sg->aBeg && sg->aBeg < p->refend){
                *prv = *fst = p;
            }else if(((pp != NULL && p == pp->next) || (pp == NULL && p == blockhead))
                        && p->refbeg > sg->aBeg && p->refbeg < sg->aEnd){
                *prv = pp;
                *fst = p;
            }else if((p->next != NULL && strcmp(p->next->refchrom,sg->aChr)==0
                        && p->refend <= sg->aBeg && sg->aBeg < p->next->refbeg)){
                *prv = p;
                *fst = p->next;
            }
            // determine lst and nxt
            if(p->refbeg < sg->aEnd && sg->aEnd <= p->refend){
                *lst = *nxt = p;
                break;
            }else if ((p->next != NULL && strcmp(p->next->refchrom,sg->aChr)==0
                        && p->refend < sg->aEnd && sg->aEnd <= p->next->refbeg)
                    || ((p->next == NULL || strcmp(p->next->refchrom,sg->aChr)!=0)
                        && p->refend < sg->aEnd)){
                *lst = p;
                *nxt = p->next;                          
                break;
            }
        }
    }
}

static struct chain_list *read_chain(string chainfile) {
	FILE *fp;
	char buf[PATH_MAX];
	struct chain_list *chainlist, *clast, *cp;
	struct gf_list *gflast, *gp;
	fp = fopen(chainfile.c_str(), "r");
	chainlist = clast = NULL;
    
	while(fgets(buf, PATH_MAX, fp)) {
		if (buf[0] == '\n' || buf[0] == '#')
			continue;
		if (buf[0] == 'c') {
			cp = (struct chain_list *)malloc(sizeof(struct chain_list));
			cp->next = NULL;	
			cp->gf = NULL;
            char fchrom[500];
            char schrom[500];
			if (sscanf(buf, "chain %*d %s %*d %c %d %d %s %d %c %d %d %d", fchrom, &(cp->forient), &(cp->fbeg), &(cp->fend), schrom, &(cp->slen), &(cp->sorient), &(cp->sbeg), &(cp->send), &(cp->cid)) != 10)
				fprintf(stderr,"cannot parse: %s", buf);
            cp->fchr = (char*)malloc(sizeof(char)*string(fchrom).length()+1);
            strcpy(cp->fchr, fchrom);
            cp->schr = (char*)malloc(sizeof(char)*string(schrom).length()+1);
            strcpy(cp->schr, schrom);
			if (chainlist == NULL)
				chainlist = clast = cp;
			else {
				clast->next = cp;
				clast = cp;
			}
		}
		else {
			gp = (struct gf_list *)malloc(sizeof(struct gf_list));
			gp->next = NULL;
			gp->fgap = gp->sgap = 0;
			if (sscanf(buf, "%d %d %d", &(gp->size), &(gp->fgap), &(gp->sgap)) == 3
				|| sscanf(buf, "%d", &(gp->size)) == 1) {
				if (cp->gf == NULL) 
					cp->gf = gflast = gp;
				else {
					gflast->next = gp;
					gflast = gp;
				}
			}else{
				fprintf(stderr,"cannot parse: %s", buf);
            }
		}
	}
	fclose(fp);
	return chainlist;

}

static void free_gf_list(struct gf_list *gflist) {
	struct gf_list *p, *q;
	p = gflist;
	for (;;) {
		q = p->next;
		free(p);
		if (q == NULL)
			break;
		else
			p = q;
	}
}

static void free_chain_list(struct chain_list *chainlist) {
	struct chain_list *p, *q;
	p = chainlist;
	for (;;) {
		q = p->next;
		free_gf_list(p->gf);
		free(p->fchr);
		free(p->schr);
		free(p);
		if (q == NULL)
			break;
		else
			p = q;
	}
}

void mapbase(int cid, string rspe, char *rchr, int rpos, 
						 string sspe, char *schr, char orient, string side,
						 int *spos, int *newrpos) {
	string chainfile;
	struct chain_list *chain ;
	struct gf_list *gfa;
	int rs, ss, ingap, roff, soff, ref, i;
    
    rs = ref_idx;
	ss = get_spc_idx(sspe);

	chainfile = chDir+"/"+spc_list[rs]+"/"+ spc_list[ss]+ "/chain/"+rchr+".chain";
	if (spctag_list[ss] == 2) {
		for (i = 0; i < spc_cnt; i++) {
			if (i == rs || i == ss){
				continue;
            }
            if (hs_exist[i]!="null"){
                free_chain_list(chain_list[i]);
                chain_list.erase(i);
                hs_exist[i]="null";
            }
		}
	}
	if (hs_exist[ss] != string(rchr)) {
		if (hs_exist[ss]!="null" ){
            free_chain_list(chain_list[ss]);
            chain_list.erase(ss);
        }
		hs_exist[ss] = string(rchr);
		chain_list[ss] = read_chain(chainfile);
	}	
	for (chain = chain_list[ss]; chain != NULL; chain = chain->next)
		if (chain->cid == cid)
			break;
	
	if (chain == NULL){
        fprintf(stderr, "[ERROR] Null chain\n");
    }
	if (rpos < chain->fbeg || rpos > chain->fend){
		fprintf(stderr, "[ERROR] Wrong reference position\n");
    }
	
	roff = soff = ingap = 0;	
	ref = rpos - chain->fbeg; //offset of ref
	
	for (gfa = chain->gf; gfa != NULL; gfa = gfa->next) {
		if (roff + gfa->size > ref) 
			break;
		else {	
			roff += gfa->size;
			soff += gfa->size;
		}
		if (roff + gfa->fgap >= ref) {
			ingap = 1;
			break;
		}
		else {
			roff += gfa->fgap;
			soff += gfa->sgap;
		}
	}
	if (ingap == 1) {
		if (side == "right") {
			roff += gfa->fgap;
			soff += gfa->sgap;
		}
	}
	else {
		soff += (ref - roff);
		roff = ref;
	}
	if (orient == '+') {
		*spos = chain->sbeg + soff;
		*newrpos = chain->fbeg + roff;
	}
	else {
		*spos = chain->slen - (chain->sbeg + soff); //rev_comp coordinate
		*newrpos = chain->fbeg + roff;
	}
}

void break_segment_position(struct seg *sg, int pos, int idx){
    struct seg *newsg;
    int lnewpos, rnewpos;
    newsg = new seg;
    newsg->next = sg->next;
    sg->next = newsg;
    newsg->aChr = (char*)malloc(sizeof(char)*string(sg->aChr).length()+1);
    strcpy(newsg->aChr,sg->aChr);
    newsg->bChr = (char*)malloc(sizeof(char)*string(sg->bChr).length()+1);
    strcpy(newsg->bChr,sg->bChr);
    newsg->aEnd = sg->aEnd;
    newsg->dir = sg->dir;
    newsg->chainID = sg->chainID;

    if (sg->dir == '+') {
        newsg->bEnd = sg->bEnd;
        mapbase(sg->chainID, ref_spc, sg->aChr, pos, spc_list[idx], sg->bChr, sg->dir, "left", &(sg->bEnd), &(lnewpos));
        mapbase(sg->chainID, ref_spc, sg->aChr, pos, spc_list[idx], sg->bChr, sg->dir, "right", &(newsg->bBeg), &(rnewpos));
    }
    else {
        newsg->bBeg = sg->bBeg;
        mapbase(sg->chainID, ref_spc, sg->aChr, pos, spc_list[idx], sg->bChr, sg->dir, "left", &(sg->bBeg), &(lnewpos));
        mapbase(sg->chainID, ref_spc, sg->aChr, pos, spc_list[idx], sg->bChr, sg->dir, "right", &(newsg->bEnd), &(rnewpos));
    }
    sg->aEnd = lnewpos;
    newsg->aBeg = rnewpos;
}
void break_block_position(struct block *blk, int pos){
    int i;
    struct block *newblk;
    struct seg *sg;

    newblk = my_allocate_newblock();
    newblk->refchrom = (char*)malloc(sizeof(char)*string(blk->refchrom).length()+1);
    strcpy(newblk->refchrom,blk->refchrom);
    newblk->refbeg = pos;
    newblk->refend = blk->refend;
    blk->refend = pos;

    for (i = 0; i < spc_cnt; i++) {
        if (i == ref_idx)
            continue;
        sg = blk->speseg[i];
        if (sg == NULL) // NULL CHECK 
            continue;
        if (pos <= sg->aBeg) {
            fill_block(newblk, i, sg);
            blk->speseg[i] = NULL;
        }
        else if (pos >= sg->aEnd) {
            continue;
        }
        else {
            break_segment_position(sg, pos, i);
            fill_block(newblk, i, sg->next);
            sg = sg->next;
            blk->speseg[i]->next = NULL;
        }
    }
    newblk->next = blk->next;
    blk->next = newblk;    
}
void add_descendant_segs(int idx, struct block **blkhead ,struct seg **sglist){
    struct seg *sg;
    struct block *blocklist, *lastblock, *newblock, *last;
    struct block *prv, *nxt, *fst, *lst;
    string prevchr = "init";
    string prevtarchr = "init";
    int cnt_chr = 0;    
    int cnt_order = 0;    
    int pos;
    blocklist = lastblock = NULL;

    if (*blkhead == NULL) {
        for (sg = *sglist; sg != NULL; sg = sg->next) {
            cnt_order ++;
            if(prevchr != string(sg->aChr)){
                prevchr = string(sg->aChr);
            }
            if(prevtarchr != string(sg->bChr)){
                prevtarchr = string(sg->bChr);
                cnt_chr++;
            }
            newblock = my_allocate_newblock();
            fill_block(newblock, idx, sg);
            if (blocklist == NULL)
                blocklist = lastblock = newblock;
            else {
                lastblock->next = newblock;
                lastblock = newblock;
            }
        }
        *blkhead = blocklist;
        maxorder = cnt_order;
    }else{
        pos = 0;
        last = NULL;
        prevchr = "init";
        prevtarchr = "init";
        cnt_chr=0;
        for(sg = *sglist; sg!=NULL ; sg = sg->next){
            cnt_order ++;
            if(prevchr != string(sg->aChr)){
                prevchr = string(sg->aChr);
            }
            if(prevtarchr != string(sg->bChr)){
                prevtarchr = string(sg->bChr);
                cnt_chr++;
            }
            find_insert_position(sg, *blkhead, last, &prv, &nxt, &fst, &lst);
            last = prv;
            if (fst == NULL && lst == NULL) {
                newblock = my_allocate_newblock();
                fill_block(newblock, idx, sg);
                if (prv == NULL && nxt == NULL) continue;
                if (prv != NULL && nxt != NULL && nxt == prv->next) {
                    newblock->next = prv->next;
                    prv->next = newblock;
                }else if(prv == NULL){
                    newblock->next = *blkhead;
                    *blkhead = newblock;
                }else if(nxt == NULL){
                    prv->next = newblock; 
                }
            }else if(fst == lst){
                if (fst->speseg[idx] == NULL){
                    fill_block(fst, idx, sg);
                }else{
                    pos = (fst->speseg[idx]->aEnd + sg->aBeg) / 2;
                    break_block_position(fst, pos);
                    fst = fst->next;
                    fill_block(fst, idx, sg);
                }
            }else{
                if (fst == NULL || lst == NULL) {
                    fprintf(stderr,"DIE2: NO INS POS:\n");
                }
                if (fst->speseg[idx] != NULL) {
                    pos = (fst->speseg[idx]->aEnd + sg->aBeg) / 2;
                    break_block_position(fst, pos);
                    fst = fst->next;
                }
                for(; fst != lst; fst = fst->next) {
                    pos = (fst->refend + fst->next->refbeg) / 2;
                    if (pos <= sg->aBeg)
                        continue;
                    break_segment_position(sg, pos, idx);
                    fill_block(fst, idx, sg);
                    sg = sg->next;
                }
                fill_block(fst, idx, sg);         
            }

        }
    }
    maxchr = maxchr < cnt_chr ? cnt_chr : maxchr;
    maxorder += maxchr;
    maxorder = maxorder < cnt_order ? cnt_order : maxorder;
}
struct block_list * PART_GENM(struct block *commonblocklist){
    int ss, rs, cflag;
    struct seg *sg;
    struct block *blk;
    struct block_list *commonblocklist4 = NULL;
    struct block_list *last_blk = new block_list();
    rs = ref_idx;
    //sanity check
    for (blk = commonblocklist; blk->next != NULL; blk = blk->next) {
        if (blk->refbeg >= blk->refend)
            fprintf(stderr,"[ERROR]end >= beg:%s.%s:%d-%d\n", ref_spc.c_str(), blk->refchrom, blk->refbeg, blk->refend);
        if (blk->refchrom == blk->next->refchrom) {
            if (blk->refend > blk->next->refbeg) {
                fprintf(stderr,"[ERROR]out of order:\n%s.%s:%d-%d %s.%s:%d-%d",ref_spc.c_str(), blk->refchrom, blk->refbeg, blk->refend, ref_spc.c_str(), blk->next->refchrom, blk->next->refbeg, blk->next->refend);
            }
        }
    }
    //print building blocks
    string build_blcks = outDir + "/Building.Blocks";
    FILE * fp3 = fopen(build_blcks.c_str(), "w");
    for (blk = commonblocklist; blk != NULL; blk = blk->next) {
        if (blk->refend - blk->refbeg <= min_rsl) continue;
        cflag = 0;
        for (ss = 0; ss < spc_cnt; ss++) {
            if (rs == ss || spctag_list[ss] != 1)
                continue;

            for (sg = blk->speseg[ss]; sg != NULL; sg = sg->next) {
                if (sg->bEnd - sg->bBeg <= min_rsl) {
                    cflag = 1;
                    break;
                }
            }
        }
        if (cflag == 1) continue;
        struct block_list *new_block;
        new_block = my_allocate_newblock_list();
        struct seg_list *ref_seg = new seg_list();
        fprintf(fp3,">\n");
        fprintf(fp3,"%s.%s:%d-%d +\n",spc_list[ref_idx].c_str(), blk->refchrom, blk->refbeg, blk->refend);
        ref_seg->chr = (char*)malloc(sizeof(char)*string(blk->refchrom).length());
        strcpy(ref_seg->chr,blk->refchrom);
        ref_seg->beg = blk->refbeg;
        ref_seg->end = blk->refend;
        ref_seg->chid = 0;
        ref_seg->orient = '+';
        new_block->speseg[rs] = ref_seg;
        for (ss = 0; ss < spc_cnt; ss++) {
            struct seg_list *p_seg = NULL;
            struct seg_list *last = new seg_list();
            if (rs == ss)
                continue;
            for (sg = blk->speseg[ss]; sg != NULL; sg = sg->next){
                fprintf(fp3, "%s.%s:%d-%d %c (%d)\n",spc_list[ss].c_str(), sg->bChr, sg->bBeg, sg->bEnd, sg->dir, sg->chainID);
                struct seg_list *p_sub = new seg_list();
                p_sub->chr = (char*)malloc(sizeof(char)*string(sg->bChr).length()+1);
                strcpy(p_sub->chr,sg->bChr);
                p_sub->beg = sg->bBeg;
                p_sub->end = sg->bEnd;
                p_sub->orient = sg->dir;
                p_sub->chid = sg->chainID;
                p_sub->next = NULL;
                if (p_seg == NULL){
                    p_seg = last = p_sub;
                }else{
                    last->next = p_sub;
                    last = p_sub;
                }
            }
            new_block->speseg[ss] = p_seg;
        }
        if(commonblocklist4 == NULL){
            commonblocklist4 = last_blk = new_block;
        }else{
            last_blk->next = new_block;
            last_blk = new_block;
        }
        fprintf(fp3,"\n");
    }
    fclose(fp3);

    return commonblocklist4;
}
void fill_block_out(struct block  *blck, int idx, struct seg *sg){
    struct seg *newsg, *p;
    newsg = new seg();
    newsg->next = NULL;
    newsg->aChr = (char*)malloc(sizeof(char)*string(sg->aChr).length()+1);
    newsg->bChr = (char*)malloc(sizeof(char)*string(sg->bChr).length()+1);
    strcpy(newsg->aChr, sg->aChr);
    strcpy(newsg->bChr, sg->bChr);
    newsg->aBeg = sg->aBeg;
    newsg->aEnd = sg->aEnd;
    newsg->bBeg = sg->bBeg;
    newsg->bEnd = sg->bEnd;
    newsg->dir = sg->dir;
    newsg->chainID = sg->chainID;

    if (blck->speseg[idx] == NULL)
        blck->speseg[idx] = newsg;
    else {
        for (p = blck->speseg[idx]; p->next != NULL; p = p->next)
            ;
        p->next = newsg;
    }
}
void add_outgroup_segs(int idx, struct block **head, struct seg **sglist) {
    struct seg *sg;
    struct block *prv, *nxt, *fst, *lst, *last;
    int pos;
    string prevchr;
    int cnt_chr = 0;    

    last = NULL;
    prevchr = "init";

    for (sg = *sglist; sg != NULL; sg = sg->next){
        if (prevchr != string(sg->aChr)) {
            prevchr = string(sg->aChr);
            cnt_chr ++;
        }
        find_insert_position(sg, *head, last, &prv, &nxt, &fst, &lst);
        last = prv;
        if (fst == NULL && lst == NULL)
            continue;
        if (fst == lst)
            fill_block_out(fst, idx, sg);
        else {
            if (fst == NULL || lst == NULL) {
                fprintf(stderr ,"DIE3: NO INS POS: \n");
            }
            for (; fst != lst; fst = fst->next) {
                pos = (fst->refend + fst->next->refbeg) / 2;
                if (pos <= sg->aBeg)
                    continue;
                break_segment_position(sg, pos, idx);
                fill_block_out(fst, idx, sg);
                sg = sg->next;
            }
            fill_block_out(fst, idx, sg);
        }
    }
    fprintf(stderr, "\n");
    maxchr = maxchr < cnt_chr ? (cnt_chr+2) : maxchr;
}
void assign_states (struct block_list *head){
    struct block_list *blk;
    struct seg_list *sg;
    int i;

    for (blk = head; blk != NULL; blk = blk->next) {
        for (i = 0; i < spc_cnt; i++) {
            if (blk->speseg[i] != NULL) {
                blk->speseg[i]->state = FIRST;
                for (sg = blk->speseg[i]->next; sg != NULL && sg->next != NULL; sg = sg->next)
                    sg->state = MIDDLE;
                if (sg != NULL)
                    sg->state = LAST;
                if (blk->speseg[i]->next == NULL)
                    blk->speseg[i]->state = BOTH;
            }
        }
    }
}
void assign_orders(struct block_list *head){
    struct block_list *blk;
    struct seg_list *sg;
    int id, subid, i;
    for (id = 0, blk = head; blk != NULL; blk = blk->next) {
        blk->id = ++id;
        for (i = 0; i < spc_cnt; i++) {
            subid = 0;
            for (sg = blk->speseg[i]; sg != NULL; sg = sg->next) {
                sg->id = blk->id;
                sg->subid = ++subid;
            }
        }
    }
}
int overlap(struct seg_list *x, struct seg_list *y) {
    int ovlp;
    int b1, e1, b2, e2, len1, len2;
    b1 = x->beg;
    e1 = x->end;
    b2 = y->beg;
    e2 = y->end;
    len1 = e1 - b1;
    len2 = e2 - b2;
    if (strcmp(x->chr,y->chr)==0 &&
            ((b1 >= b2 && e1 <= e2) || (b1 <= b2 && e1 >= e2)
                || (b1 < b2 && e1 > b2 && e1 - b2 > MINOVL * min(len1, len2))
                || (b1 < e2 && e1 > e2 && e2 - b1 > MINOVL * min(len1, len2)) ))
        ovlp = 1;
    else
        ovlp = 0;
    return ovlp;
}

int contain_all(struct block_list *blst){
    int i;
    for (i = 0; i < spc_cnt; i++)
        if (spctag_list[i] == 1 && blst->speseg[i] == NULL)
            break;
    return (i == spc_cnt) ? 1 : 0;

}
void clean_up(struct block_list **head){
    struct block_list *p, *q;
    int i, lenp, lenq, rs;
    rs = ref_idx;
    for (p = *head; p != NULL; p = p->next) {
        if (!contain_all(p))
            continue;
        for (q = p->next; q != NULL; q = q->next) {
            if (!contain_all(q))
                continue;
            for (i = 0; i < spc_cnt; i++)
                if (spctag_list[i] == 1 && !overlap(p->speseg[i], q->speseg[i]))
                    break;
            lenp = p->speseg[rs]->end - p->speseg[rs]->beg;
            lenq = q->speseg[rs]->end - q->speseg[rs]->beg;
            if (i == spc_cnt) {
                if (lenp < lenq)
                    p->isdup = 1;
                else
                    q->isdup = 1;
            }
        }
    }

    q = NULL;
    for (p = *head; p != NULL;) {
        if (p->isdup) {
            if (q == NULL) {
                q = p;
                p = p->next;
                *head = p;
                q->next = NULL;
                free_block_list(q);
                q = NULL;
            }
            else {
                q->next = p->next;
                p->next = NULL;
                free_block_list(p);
                p = q->next;
            }
        }
        else {
            q = p;
            p = p->next;
        }
    }

}
int with_longblock(struct block_list *blk){
    struct block_list *p;
    int i, rslen, tarlen, illegal;

    illegal = 0;
    p = blk;
    rslen = p->speseg[ref_idx]->end - p->speseg[ref_idx]->beg;

    for (i = 0; i < spc_cnt; i++) {
        if (spctag_list[i] == 1 && p->speseg[i] != NULL) {
            tarlen = p->speseg[i]->end - p->speseg[i]->beg;
            if (tarlen > rslen * 2) break;
        }
    }
    if (i != spc_cnt)
        illegal = 1;

    return illegal;

}
void clean_up_long(struct block_list **blockhead){
    struct block_list *p, *q;
    q = NULL;
    for (p = *blockhead; p != NULL;) {
        if (with_longblock(p)) {
            if (q == NULL) {
                q = p;
                p = p->next;
                *blockhead = p;
                q->next = NULL;
                free_block_list(q);
                q = NULL;
            }
            else {
                q->next = p->next;
                p->next = NULL;
                free_block_list(p);
                p = q->next;
            }
        }
        else {
            q = p;
            p = p->next;
        }
    }
}
int random_piece(struct seg_list *sg) {
    int ru = 0;
    if (strstr(sg->chr, "chrUn") != NULL || strstr(sg->chr, "random") || strstr(sg->chr, "chrY") != NULL || strstr(sg->chr, "chrM") != NULL)
        ru = 1;
    return ru;
}
int messy_piece(struct seg_list *sg, struct block_list *blk, int idx) {
    struct block_list *b;
    struct seg_list *p;
    int b1, e1, b2, e2, len1, len2, messy;
    b1 = sg->beg;
    e1 = sg->end;
    len1 = e1 - b1;
    messy = 0;
    for (b = blk; b != NULL; b = b->next) {
        if ((p = b->speseg[idx]) != NULL ) {
            if (p == sg)
                continue;
            b2 = p->beg;
            e2 = p->end;
            len2 = e2 - b2;
            if (strcmp(p->chr,sg->chr)==0 &&
                    ((b1 >= b2 && e1 <= e2)
                || (b1 <= b2 && e1 <= e2 && e1 > b2 && b2 - b1 < AFEW * len1 && len1 <= len2)
                || (b1 >= b2 && e1 >= e2 && b1 < e2 && e1 - e2 < AFEW * len1 && len1 <= len2))) {
                messy = 1;
                break;
            }
        }
    }
    return messy;
}
void clean_up_again(struct block_list *head){
    struct block_list *p;
    struct seg_list *sg, *tg;
    int i;
    for (p = head; p != NULL; p = p->next) {
        for (i = 0; i < spc_cnt; i++) {
            if (i == ref_idx){
                continue;
            }
            for (sg = p->speseg[i]; sg != NULL; ) {
                    if (random_piece(sg) || messy_piece(sg, head, i)) {
                        if (sg == p->speseg[i]) {
                            p->speseg[i] = sg->next;
                            tg = sg = p->speseg[i];
                        }else {
                            tg->next = sg->next;
                            sg = tg->next;
                        }
                    }else {
                        tg = sg;
                        sg = sg->next;
                }
            }//for
        }//for
    }//for
}
int illegal_block(struct block_list *blk) {
    struct block_list *p;
    int i, len, illegal;

    illegal = 0;
    p = blk;
    len = p->speseg[ref_idx]->end - p->speseg[ref_idx]->beg;
    if (len < min_rsl) {
        illegal = 1;
    }
    for (i = 0; i < spc_cnt; i++) {
        if (spctag_list[i] == 1 && (p->speseg[i] == NULL
                || p->speseg[i]->end - p->speseg[i]->beg < len * MINDESSEG)) {
            break;
        }
    }
    if (i != spc_cnt)
        illegal = 1;

    return illegal;
}
void trim(struct block_list **blockhead){
    struct block_list *p, *q;
    q = NULL;
    for (p = *blockhead; p != NULL;) {
        if (illegal_block(p)) {
            if (q == NULL) {
                q = p;
                p = p->next;
                *blockhead = p;
                q->next = NULL;
                free_block_list(q);
                q = NULL;
            }
            else {
                q->next = p->next;
                p->next = NULL;
                free_block_list(p);
                p = q->next;
            }
        }
        else {
            q = p;
            p = p->next;
        }
    }
}
struct block_list * MKOR_BLKS(struct block_list **commonblocklist4){
    struct block_list *commonblocklist = *commonblocklist4;
    assign_states(commonblocklist);
    assign_orders(commonblocklist);
    clean_up(&commonblocklist);
    clean_up_long(&commonblocklist);
    clean_up_again(commonblocklist);
    trim(&commonblocklist);
    assign_states(commonblocklist);
    assign_orders(commonblocklist);
    string orthology_blcks = outDir + "/_orthology.blocks.tmp";
    FILE * fp4 = fopen(orthology_blcks.c_str(), "w");
    
    int i;
    struct seg_list *sg;
    struct block_list *blist;
    for (blist = commonblocklist; blist != NULL; blist = blist->next) {
        fprintf(fp4, ">%d\n", blist->id);
        for (sg = blist->speseg[ref_idx]; sg != NULL; sg = sg->next){
            fprintf(fp4, "%s.%s:%d-%d %c [%d] (%d)\n", ref_spc.c_str(), sg->chr, sg->beg, sg->end, sg->orient, sg->state, sg->chid);

        }
        for (int i = 0; i < spc_list.size(); i++) {
            if (i == ref_idx){continue;}
            string spc = spc_list[i];
            for (sg = blist->speseg[i]; sg != NULL; sg = sg->next){
                fprintf(fp4, "%s.%s:%d-%d %c [%d] (%d)\n", spc.c_str(), sg->chr, sg->beg, sg->end, sg->orient, sg->state, sg->chid);
            }
        }
        fprintf(fp4, "\n");
    }
    fclose(fp4);
    string awk_cmd = "awk '{if (NF > 2) {print $1, $2} else {print $0}}' " + orthology_blcks + "> " + outDir + "/Orthology.Blocks";
    system(awk_cmd.c_str());
	return commonblocklist;
}
void ORTH_TODR(struct block_list **commonblocklist4, int **perm){
    struct block_list *blkhead = *commonblocklist4;
    struct block_list *commonblocklist = *commonblocklist4;

    int i, j;
    struct seg_list ***head;
    struct seg_list *pp, *p, *sg, *q;
    struct block_list  *bk;
    int total[spc_cnt];
    head = (struct seg_list ***)malloc(sizeof(struct seg_list**) * spc_cnt);
    for (j = 0; j < spc_cnt; j++) {
        total[j] = 0;
        head[j] = (struct seg_list **)malloc(sizeof(struct seg_list*) * maxchr);
        for (i = 0; i < maxchr; i++){
            head[j][i] = NULL;
        }
    }
    
    for (bk = blkhead; bk != NULL; bk = bk->next) {
        for (i = 0; i < spc_cnt; i++) {
            if (spctag_list[i] == 2){
                continue;
            }
            sg = bk->speseg[i];
            for (j = 0; j < total[i]; j++){
                if (strcmp(head[i][j]->chr,sg->chr)==0){
                    break;
                }
            }
            q = new seg_list();
            q->cidlist = NULL;
            q->next = NULL;
            q->id = sg->id;
            q->beg = sg->beg;
            q->end = sg->end;
            q->orient = sg->orient;
            q->chr = (char*)malloc(sizeof(char)*string(sg->chr).length()+1);
            strcpy(q->chr, sg->chr);
            
            if (j == total[i]) {
                head[i][j] = q;
                ++total[i];
            }else{
                pp = NULL;
                for (p = head[i][j]; p != NULL; p = p->next) {
                    if ((pp == NULL && q->beg < p->beg) || (pp != NULL && q->beg < p->beg && q->beg > pp->beg)){
                        break;
                    }
                    pp = p;
                }
                if (pp == NULL) {
                    q->next = p;
                    head[i][j] = q;
                }else {
                    q->next = pp->next;
                    pp->next = q;
                }
            }

        }
    }
    string order_ds = outDir + "/_order.DS";
    FILE * fp5 = fopen(order_ds.c_str(), "w");
    for (i = 0; i < spc_cnt; i++) {
        int j = 1;
        if (spctag_list[i] == 2){
            continue;
        }
        fprintf(fp5, ">%s\n", spc_list[i].c_str());
        int k = 0 ; 
        for (j = 0; j < total[i]; j++) {
            fprintf(fp5, "# %s\n", head[i][j]->chr);
            for (p = head[i][j]; p != NULL; p = p->next) {
                if (p->orient == '+'){
                    fprintf(fp5, "%d ", p->id);
                    perm[i][k] = p->id;
                }else{
                    fprintf(fp5, "-%d ", p->id);
                    perm[i][k] = -1*(p->id);
                }
                k++;
            }
            fprintf(fp5, "$\n");
            k++;
        }
        fprintf(fp5, "\n");
    }
    fclose(fp5);
    for (j = 0; j < spc_cnt; j++) {
        for (i = 0; i < maxchr; i++){
            free_seg_list(head[j][i]);
        }
    }
}
void merge_blocks(struct block_list **blkhead, int start, int terminal){
	struct block_list *p, *q;
	struct seg_list *b;
	int i, j;
	if (terminal < start)
		fprintf(stderr, "DIE: start >terminal %d %d", start, terminal);
	for (p = *blkhead; p != NULL; p = p->next) {
        if (p->id == start){
			break;
        }
    }
	if (start == terminal) {
		for (i = 0; i < spc_cnt; i++) {
			if (spctag_list[i] == 1) {
				b = p->speseg[i];
				b->chnum = 1;
				b->cidlist = (int *)malloc(sizeof(int) * b->chnum);
				b->cidlist[0] = b->chid;
			}
		}
		return;
	}

	for (q = p; q != NULL && q->id <= terminal; q = q->next) {
		for (i = 0; i < spc_cnt; i++) {
			if (spctag_list[i] == 2)
				continue;
			if (q->speseg[i]->next != NULL)
				fprintf(stderr, "DIE: illegal block %d", q->id);
		}
	}

	for (i = 0; i < spc_cnt; i++) {
		q = p;
		if (spctag_list[i] == 1) {
			b = q->speseg[i];
			b->chnum = terminal - start + 1;
			b->cidlist = (int *)malloc(sizeof(int) * b->chnum);
			j = 0;
			while (q != NULL && q->id <= terminal) {
				b->cidlist[j++] = q->speseg[i]->chid;
				q = q->next;
			}
		}
	}

	q = p->next;
	while (q != NULL && q->id <= terminal) {
		p->next = q->next;
		q->next = NULL;
		for (i = 0; i < spc_cnt; i++) {
			if (spctag_list[i] == 2)
				continue;
			p->speseg[i]->beg = min(p->speseg[i]->beg, q->speseg[i]->beg);
			p->speseg[i]->end = max(p->speseg[i]->end, q->speseg[i]->end);
		}
		for (i = 0; i < spc_cnt; i++) {
			if (spctag_list[i] != 2)
				continue;
			if (p->speseg[i] == NULL)
				p->speseg[i] = q->speseg[i];
			else {
				for (b = p->speseg[i]; b->next != NULL; b = b->next)
					;
				b->next = q->speseg[i];
			}
		}
		for (i = 0; i < spc_cnt; i++)
			q->speseg[i] = NULL;
		free_block_list(q);
		q = p->next;
	}
}
void MKCN_SGMT(struct block_list **commonblocklist4, int **perm){
    struct block_list *blkhead = *commonblocklist4;
    struct block_list *blklast, *s;
    int i, j, rs, num, total, k, start, terminal, count;
    int status[spc_cnt];
    char *pt;
    struct seg_list *p;
    total = count = 0;
    for (s = blkhead; s != NULL; s = s->next){
        ++total;
    }

    start = terminal = 1;
    while (terminal <= total) {
        for (i = 0; i < spc_cnt; i++){
            status[i] = 0;
        }
        for (i = 0; i < spc_cnt; i++) {
            if (spctag_list[i] == 2){
                continue;
            }
            for (k = 0; k < maxorder; k++){
                if (abs(perm[i][k]) == terminal){
                    break;
                }
            }
            if ((perm[i][k] > 0 && perm[i][k+1] == terminal + 1) || (perm[i][k] < 0 && perm[i][k-1] == -terminal - 1)){
                status[i] = 1;
            }
        }
        for (i = 0; i < spc_cnt; i++) {
            if (i == ref_idx || spctag_list[i] == 2){continue;}
            status[ref_idx] = status[ref_idx] & status[i];
        }
        if (status[ref_idx] == 1){
            terminal ++;
        }else{
            merge_blocks(&blkhead, start, terminal);
            start = terminal + 1;
            terminal = start;
        }
    }
    assign_states(blkhead);
    assign_orders(blkhead);
    string conserved_segs = outDir + "/_conserved.segments.tmp";
    FILE * fp6 = fopen(conserved_segs.c_str(), "w");
    for (s = blkhead; s != NULL; s = s->next) {
		fprintf(fp6, ">%d\n", s->id);
		// for reference species first
		for (p = s->speseg[ref_idx]; p != NULL; p = p->next) {
			fprintf(fp6, "%s.%s:%d-%d %c [%d]", spc_list[ref_idx].c_str(), p->chr,p->beg, p->end, p->orient, p->state);
			if (spctag_list[ref_idx] == 0) {
				fprintf(fp6, "\n");
				continue;
			}
			else if (spctag_list[ref_idx] == 1) {
				fprintf(fp6, " {%d", p->chnum);
				for (j = 0; j < p->chnum; j++)
					printf(",%d", p->cidlist[j]);
				fprintf(fp6, "}\n");
			}
			else {
				fprintf(fp6, " (%d)\n", p->chid);
			}
		}
		
		for (i = 0; i < spc_cnt; i++) {
			if (ref_idx == i) { continue; }
			for (p = s->speseg[i]; p != NULL; p = p->next) {
				fprintf(fp6, "%s.%s:%d-%d %c [%d]", spc_list[i].c_str(), p->chr, p->beg, p->end, p->orient, p->state);
				if (spctag_list[i] == 0) {
					printf("\n");
					continue;
				}
				else if (spctag_list[i] == 1) {
					fprintf(fp6, " {%d", p->chnum);
					for (j = 0; j < p->chnum; j++){
						fprintf(fp6, ",%d", p->cidlist[j]);
                    }
					fprintf(fp6, "}\n");
				}
				else {
					fprintf(fp6, " (%d)\n", p->chid);
				}
			}
		}
		fprintf(fp6, "\n");
	}
    fclose(fp6);
}
void OUTS_ORDR(struct block_list **commonblocklist4,int **perm_id,int **perm_sid, int *outorder){
    struct block_list *blkhead = *commonblocklist4;
    int i, j;
    struct seg_list ***head;
    struct seg_list *pp, *p, *q, *sg;
    struct block_list  *bk;
    int total[spc_cnt];

    head = (struct seg_list ***)malloc(sizeof(struct seg_list**) * spc_cnt);
    for (j = 0; j < spc_cnt; j++) {
        total[j] = 0;
        head[j] = (struct seg_list **)malloc(sizeof(struct seg_list*) * maxchr);
        for (i = 0; i < maxchr; i++){
            head[j][i] = NULL;
        }
    }
    for (bk = blkhead; bk != NULL; bk = bk->next) {
        for (i = 0; i < spc_cnt; i++) {
            if (spctag_list[i] != 2){
                continue;
            }
            for(sg = bk->speseg[i]; sg!=NULL ; sg = sg->next){
                for (j = 0; j < total[i]; j++){
                    if (strcmp((head[i][j])->chr,sg->chr)==0){
                        break;
                    }
                }
                q = new seg_list();
                q->cidlist = NULL;
                q->next = NULL;
                q->id = sg->id;
                q->subid = sg->subid;
                q->beg = sg->beg;
                q->end = sg->end;
                q->orient = sg->orient;
                q->chr = (char*)malloc(sizeof(char)*string(sg->chr).length()+1);
                strcpy(q->chr, sg->chr);

                if (j == total[i]) {
                    head[i][j] = q;
                    ++total[i];
                }else{
                    pp = NULL;
                    for (p = head[i][j]; p != NULL; p = p->next) {
                        if ((pp == NULL && q->beg < p->beg) || (pp != NULL && q->beg < p->beg && q->beg > pp->beg)){
                            break;
                        }
                        pp = p;
                    }
                    if (pp == NULL) {
                        q->next = p;
                        head[i][j] = q;
                    }else {
                        q->next = pp->next;
                        pp->next = q;
                    }
                }
            }
        }
    }
    string order_og = outDir + "/_order.OG";
    int k = 0 ;
    FILE * fp7 = fopen(order_og.c_str(), "w");
    for (i = 0; i < spc_cnt; i++) {
        if (spctag_list[i] != 2){
            continue;
        }
        fprintf(fp7, ">%s\n", spc_list[i].c_str());
        for (j = 0; j < total[i]; j++) {
            fprintf(fp7, "# %s\n", head[i][j]->chr);
            for (p = head[i][j]; p != NULL; p = p->next) {
                k = ++ outorder[i];
                if (p->orient == '+'){
                    fprintf(fp7, "%d.%d ", p->id, p->subid);
                    perm_id[i][k] = p->id;
                    perm_sid[i][k] = p->subid;
                }else{
                    fprintf(fp7, "-%d.%d ", p->id, p->subid);
                    perm_id[i][k] = -1*p->id;
                    perm_sid[i][k] = p->subid;
                }
            }
            fprintf(fp7, "$\n");
            ++ outorder[i];
            ++ outorder[i];
        }
        fprintf(fp7, "\n");
    }
    fclose(fp7);
}
void merge_segs(struct block_list *blkhead, int id, int ss, int start, int terminal){
    struct block_list *p;
    struct seg_list *b, *nb;
    int j;

    if (terminal < start)
        fprintf(stderr,"DIE: start > terminal %d %d", start, terminal);

    for (p = blkhead; p != NULL; p = p->next){
        if (p->id == id){
            break;
        }
    }

    for (b = p->speseg[ss]; b != NULL; b = b->next){
        if (b->subid == start){
            break;
        }
    }

    if(b == NULL){
        return;
    }
    if(start == terminal) {
        b->chnum = 1;
        b->cidlist = (int *)malloc(sizeof(int) * b->chnum);
        b->cidlist[0] = b->chid;
        return;
    }

    b->chnum = terminal - start + 1;
    b->cidlist = (int *)malloc(sizeof(int) * b->chnum);
    nb = b;
    j = 0;
    while (nb != NULL && nb->subid <= terminal) {
        b->cidlist[j++] = nb->chid;
        nb = nb->next;
    }

    for (nb = b->next; nb != NULL; ) {
        if (nb->subid <= terminal) {
            b->beg = min(b->beg, nb->beg);
            b->end = max(b->end, nb->end);
            b->next = nb->next;
            nb->next = NULL;
            free_seg_list(nb);
            nb = b->next;
        }
        else
            break;
    }
}
void remove_tiny_pieces(struct block_list *head) {
  struct block_list *p;
    struct seg_list *sg, *tg;
    int i, len, reflen, rs;
    rs = ref_idx;
    for (p = head; p != NULL; p = p->next) {
        reflen = p->speseg[rs]->end - p->speseg[rs]->beg;
        for (i = 0; i < spc_cnt; i++) {
            if (spctag_list[i] != 2)
                continue;
            for (sg = p->speseg[i]; sg != NULL; ) {
                len = sg->end - sg->beg;
                if (len < MINOUTSEG * reflen) {
                    if (sg == p->speseg[i]) {
                        p->speseg[i] = sg->next;
                        sg->next = NULL;
                        free_seg_list(sg);
                        tg = sg = p->speseg[i];
                    }
                    else {
                        tg->next = sg->next;
                        sg->next = NULL;
                        free_seg_list(sg);
                        sg = tg->next;
                    }
                }
                else {
                    tg = sg;
                    sg = sg->next;
                }
            }
        }
    }
}
void merge_chlist(struct block_list *head) {
    struct block_list *blk;
    struct seg_list *sg;
    int buf[5000], j, prev, i, k;

    for (blk = head; blk != NULL; blk = blk->next) {
        for (i = 0; i < spc_cnt; i++) {
            if (spctag_list[i] == 0)
                continue;
            for (sg = blk->speseg[i]; sg != NULL; sg = sg->next) {
                prev = j = 0;
                for (k = 0; k < sg->chnum; k++) {
                    if (sg->cidlist[k] != prev) {
                        buf[j++] = sg->cidlist[k];
                        prev = sg->cidlist[k];
                    }
                }
                if (j != sg->chnum) {
                    sg->chnum = j;
                    free(sg->cidlist);
                    sg->cidlist = (int *)malloc(sizeof(int) * j);
                    for (k = 0; k < j; k++)
                        sg->cidlist[k] = buf[k];
                }
            }
        }
    }

}
void CLEN_OGSG(struct block_list **commonblocklist4,int **perm_id,int **perm_sid, int *outorder){
    struct block_list *blkhead = *commonblocklist4;
    struct block_list *bk;
    int i, j, rs, num, snum, total, k, start, terminal;
    char *pt;
    struct seg_list *p;
    for (total = 0, bk = blkhead; bk != NULL; bk = bk->next){
        ++total;
    }
    for (i = 0 ; i < spc_cnt; i ++){
        if (spctag_list[i] != 2){
            continue;
        }
        start = terminal = 1;
        for (j = 1 ; j <= total ; j ++ ){
            for (;;){
                for ( k = 0 ; k <= outorder[i] ; k ++){
                    if (abs(perm_id[i][k]) == j && perm_sid[i][k] == terminal){
                        break;
                    }
                }
                if (k > outorder[i]){
                    start = terminal = 1;
                    break;
                }
                if((perm_id[i][k] > 0 && perm_id[i][k+1] == perm_id[i][k] && perm_sid[i][k+1] == terminal + 1)
                    || (perm_id[i][k] < 0 && perm_id[i][k-1] == perm_id[i][k] && perm_sid[i][k-1] == terminal + 1)){
                    terminal ++;
                }else{
                    merge_segs(blkhead, j, i, start, terminal);
                    start = terminal + 1;
                    terminal = start;

                }
            }
        }
    }
    remove_tiny_pieces(blkhead);
    assign_states(blkhead);
    merge_chlist(blkhead);
    string conserved_segs2 = outDir + "/_conserved.segments.tmp2";
    FILE * fp8 = fopen(conserved_segs2.c_str(), "w");
    for (bk = blkhead; bk != NULL; bk = bk->next) {
        fprintf(fp8, ">%d\n", bk->id);
        // for reference species first
        for (p = bk->speseg[rs]; p != NULL; p = p->next) {
            fprintf(fp8, "%s.%s:%d-%d %c [%d] [%d.%d]", spc_list[ref_idx].c_str(), p->chr, p->beg, p->end, p->orient, p->state, p->id, p->subid);
            if (spctag_list[ref_idx] == 0){
                fprintf(fp8, "\n");
                continue;
            }
            fprintf(fp8, " {%d", p->chnum);
            for (j = 0; j < p->chnum; j++){
                fprintf(fp8, ",%d", p->cidlist[j]);
            }
            fprintf(fp8, "}\n");
        }

        for (i = 0; i < spc_cnt; i++) {
            if (ref_idx == i) { continue; }
            for (p = bk->speseg[i]; p != NULL; p = p->next) {
                fprintf(fp8, "%s.%s:%d-%d %c [%d] [%d.%d]", spc_list[i].c_str(), p->chr, p->beg, p->end, p->orient, p->state, p->id, p->subid);
                if (spctag_list[i] == 0){
                    fprintf(fp8, "\n");
                    continue;
                }
                fprintf(fp8, " {%d", p->chnum);
                for (j = 0; j < p->chnum; j++){
                    fprintf(fp8, ",%d", p->cidlist[j]);
                }
                fprintf(fp8, "}\n");
            }
        }
        fprintf(fp8, "\n");
    }
    fclose(fp8);
    string awk_cmd = "awk '{if (NF > 2) {print $1, $2} else {print $0}}' " + conserved_segs2 + "> " + outDir + "/Conserved.Segments";
    system(awk_cmd.c_str());
}
void print_join(FILE *fp, int l, char ol, int r, char oor) {
    if (l == r)
        return;
    if (ol == '-')
        l = -l;
    fprintf(fp, "%*d", 5, l);
    if (oor == '-')
        r = -r;
    fprintf(fp, "\t%*d\n", 5, r);
}
void CRET_GNOM(struct block_list **commonblocklist4){
    FILE *jfp;
    int i, j, count, rs;
    char buf[500];
    struct seg_list ***head;
    struct seg_list *pp, *p, *q, *sg;
    struct block_list *blkhead, *bk;
    int total[spc_cnt];
    head = (struct seg_list ***)malloc(sizeof(struct seg_list**) * spc_cnt);
    for (j = 0; j < spc_cnt; j++) {
        total[j] = 0;
        head[j] = (struct seg_list **)malloc(sizeof(struct seg_list*) * maxchr);
        for (i = 0; i < maxchr; i++){
            head[j][i] = NULL;
        }
    }
    j = 0;
    blkhead = *commonblocklist4;
    for (bk = blkhead; bk != NULL; bk = bk->next)  {
        for (i = 0; i < spc_cnt; i++) {
            for (sg = bk->speseg[i]; sg != NULL; sg = sg->next) {
                for (j = 0; j < total[i]; j++){
                    if (strcmp((head[i][j])->chr,sg->chr)==0){
                        break;
                    }
                }
                q = new seg_list;
                q->cidlist = NULL;
                q->next = NULL;
                q->id = sg->id;
                q->subid = sg->subid;
                q->beg = sg->beg;
                q->end = sg->end;
                q->orient = sg->orient;
                q->state = sg->state;
                q->chr = (char*)malloc(sizeof(char)*string(sg->chr).length()+1);
                strcpy(q->chr, sg->chr);

                if (j == total[i]) {
                    head[i][j] = q;
                    ++total[i];
                }else {
                    pp = NULL;
                    for (p = head[i][j]; p != NULL; p = p->next) {
                        if ((pp == NULL && q->beg < p->beg) || (pp != NULL && q->beg < p->beg && q->beg > pp->beg))
                            break;
                        pp = p;
                    }
                    if (pp == NULL) {
                        q->next = p;
                        head[i][j] = q;
                    }
                    else {
                        q->next = pp->next;
                        pp->next = q;
                    }
                }
            }
        }
    }
    for (count = 0, bk = blkhead; bk != NULL; bk = bk->next){
		++count;
    }
    // for reference species first
    string order_GO = outDir + "/Genomes.Order";
    FILE * fp9 = fopen(order_GO.c_str(), "w");
	fprintf(fp9, ">%s\t%d\n", spc_list[ref_idx].c_str(), total[ref_idx]);
	for (j = 0; j < total[ref_idx]; j++) {
		fprintf(fp9, "# %s\n", (head[ref_idx][j])->chr);
		for (p = head[ref_idx][j]; p != NULL; p = p->next) {
			if (p->orient == '+')
				fprintf(fp9,"%d ", p->id);
			else
				fprintf(fp9,"-%d ", p->id);
		}
		fprintf(fp9,"$\n");
	}
	fprintf(fp9,"\n");

	for (i = 0; i < spc_cnt; i++) {
		if (spctag_list[i] == 2 || ref_idx == i)
			continue;
		fprintf(fp9,">%s\t%d\n", spc_list[i].c_str(), total[i]);
		for (j = 0; j < total[i]; j++) {
			fprintf(fp9,"# %s\n", (head[i][j])->chr);
			for (p = head[i][j]; p != NULL; p = p->next) {
				if (p->orient == '+')
					fprintf(fp9,"%d ", p->id);
				else
					fprintf(fp9,"-%d ", p->id);
			}
			fprintf(fp9,"$\n");
		}
		fprintf(fp9,"\n");
	}
    fclose(fp9);

	// for reference species first
	string ref_jo = outDir+"/"+spc_list[ref_idx]+".joins";
	jfp = fopen(ref_jo.c_str(), "w");
	fprintf(jfp, "#%d\n", count);
	for (j = 0; j < total[ref_idx]; j++) {
		p = head[ref_idx][j];
		if (p == NULL)
			continue;
		if (spctag_list[ref_idx] != 2)
			if ((p->state == FIRST && p->orient == '+') || p->state == BOTH
					|| (p->state == LAST && p->orient == '-'))
				print_join(jfp, 0, '+', p->id, p->orient);
		for (; p->next != NULL; p = p->next) {
			q = p->next;
			if (((p->state == FIRST && p->orient == '-') || (p->state == LAST && p->orient == '+'))
					&& ((q->state == FIRST && q->orient == '+') || (q->state == LAST && q->orient == '-')))
				print_join(jfp, p->id, p->orient, q->id, q->orient);
			if (p->state == BOTH && ((q->state == FIRST && q->orient == '+')
				||(q->state == LAST && q->orient == '-')))
				print_join(jfp, p->id, p->orient, q->id, q->orient);
			if (((p->state == FIRST && p->orient == '-') || (p->state == LAST && p->orient == '+'))
					&& q->state == BOTH)
				print_join(jfp, p->id, p->orient, q->id, q->orient);
			if (p->state == BOTH && q->state == BOTH)
				print_join(jfp, p->id, p->orient, q->id, q->orient);
		}
		if (spctag_list[ref_idx] != 2)
			if (p->state == BOTH || (p->state == LAST && p->orient == '+')
					|| (p->state == FIRST && p->orient == '-'))
				print_join(jfp, p->id, p->orient, 0, '+');
	}
	fclose(jfp);

	for (i = 0; i < spc_cnt; i++) {
		if (ref_idx == i) { continue; }
        ref_jo = outDir+"/"+spc_list[i]+".joins";
	    jfp = fopen(ref_jo.c_str(), "w");
	    fprintf(jfp, "#%d\n", count);

		for (j = 0; j < total[i]; j++) {
			p = head[i][j];
			if (p == NULL)
				continue;
			if (spctag_list[i] != 2)
				if ((p->state == FIRST && p->orient == '+') || p->state == BOTH
						|| (p->state == LAST && p->orient == '-'))
					print_join(jfp, 0, '+', p->id, p->orient);
			for (; p->next != NULL; p = p->next) {
				q = p->next;
				if (((p->state == FIRST && p->orient == '-') || (p->state == LAST && p->orient == '+'))
						&& ((q->state == FIRST && q->orient == '+') || (q->state == LAST && q->orient == '-')))
					print_join(jfp, p->id, p->orient, q->id, q->orient);
				if (p->state == BOTH && ((q->state == FIRST && q->orient == '+')
					||(q->state == LAST && q->orient == '-')))
					print_join(jfp, p->id, p->orient, q->id, q->orient);
				if (((p->state == FIRST && p->orient == '-') || (p->state == LAST && p->orient == '+'))
						&& q->state == BOTH)
					print_join(jfp, p->id, p->orient, q->id, q->orient);
				if (p->state == BOTH && q->state == BOTH)
					print_join(jfp, p->id, p->orient, q->id, q->orient);
			}
			if (spctag_list[i] != 2)
				if (p->state == BOTH || (p->state == LAST && p->orient == '+')
						|| (p->state == FIRST && p->orient == '-'))
					print_join(jfp, p->id, p->orient, 0, '+');
		}
		fclose(jfp);
	}
}
int main(int argc, char * argv[]){
    string cur_exec_name = argv[0];
    vector < string > argList(argv, argv + argc);
    string f_config = argList[1];
    char curDir[PATH_MAX];
	getcwd(curDir,PATH_MAX);
	outDir = string(curDir);

    READ_CONFIG(f_config);

    for(int i = 0 ; i < spc_cnt; i++){
        hs_exist[i]= "null";
    }
    struct block *commonblocklist = NULL;
    struct block_list *commonblocklist4 = NULL;
    //ingroup
    cerr<<"Process_chainNet " << endl;
    for (int i = 0; i < spc_list.size(); i++) {
        if (spc_list[i] == ref_spc) {
            continue;
        }
        if (spctag_list[i] == 1) {
            struct seg *seglst;
        	string tmpNetDir = netDir + "/" + ref_spc + "/" + spc_list[i] + "/net";
        	READ_NETS(i, tmpNetDir);
            seglst = PROC_SEGS(i);
            add_descendant_segs(i, &commonblocklist ,&seglst);
            free_seg(seglst);
        }else{continue;}
    }
    //outgroup
    for (int i = 0; i < spc_list.size(); i++) {
        if (spc_list[i] == ref_spc) {
            continue;
        }
        if (spctag_list[i] == 2) {
            struct seg *seglst;
        	string tmpNetDir = netDir + "/" + ref_spc + "/" + spc_list[i] + "/net";
        	READ_NETS(i, tmpNetDir);
            seglst  = PROC_SEGS(i);
            add_outgroup_segs(i, &commonblocklist ,&seglst);
            free_seg(seglst);
        }else{continue;}
    }
    cerr<<"Partition Genomes " << endl;
    commonblocklist4 = PART_GENM(commonblocklist);
    free_block(commonblocklist);
    cerr<<"Make Orthology Blocks " <<endl;
    commonblocklist4 = MKOR_BLKS(&commonblocklist4);
    int **perm = (int **)malloc(sizeof(int*) * spc_cnt);
    for (int i = 0; i < spc_cnt; i++) {
        perm[i] = (int *)malloc(sizeof(int) * maxorder);
        if (spctag_list[i] == 2){
            continue;
        }
        for (int j = 0; j < maxorder; j++){
            perm[i][j] = 0;
        }
    }

    cerr<<"Orthology Blocks To Order " << endl; 
    ORTH_TODR(&commonblocklist4, perm);
    cerr<<"Make Conseerved Segments " << endl;
    MKCN_SGMT(&commonblocklist4, perm);

    int *outorder= (int*)malloc(sizeof(int*)*spc_cnt);
    int **perm_id = (int **)malloc(sizeof(int*) * spc_cnt);
    int **perm_sid = (int **)malloc(sizeof(int*) * spc_cnt);
    for (int i = 0; i < spc_cnt; i++) {
        perm_id[i] = (int *)malloc(sizeof(int) * (maxorder*10));
        perm_sid[i] = (int *)malloc(sizeof(int) * (maxorder*10));
        outorder[i] = 0;
        if (spctag_list[i] != 2){
            continue;
        }
        for (int j = 0; j < maxorder; j++){
            perm_id[i][j] = 0;
            perm_sid[i][j] = 0;
        }
    }
    cerr<<"Order Outgroup Segments " << endl; 
    OUTS_ORDR(&commonblocklist4, perm_id,perm_sid, outorder);
    cerr<<"Clean Outgroup Segments " << endl;
    CLEN_OGSG(&commonblocklist4, perm_id, perm_sid, outorder);
    cerr<<"Create Genomes " << endl; 
    CRET_GNOM(&commonblocklist4);
    cerr<<"All Done " << endl;


    return 0;
}
