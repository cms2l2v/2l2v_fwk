#ifndef _jsonwrapper_cc_
#define _jsonwrapper_cc_

#include <string.h>
#include <string>
#include <vector>
#include <iostream>
#include <iterator>

namespace JSONWrapper{

  std::string removeWhiteSpace(const std::string& in, int N=-1);
  std::string removeQuotes(const std::string& in);
  size_t findComma(const std::string& in, int start);
  size_t findEndBrace(const std::string& in, int start);
  size_t findEndBracket(const std::string& in, int start);
  bool isObject(const std::string& in);
  bool isArray (const std::string& in);
  bool isComma (const std::string& in);

  class Object
  {
  public:
    Object(){}
    Object(const std::string& in, bool isInputFile=false, int length=-1);
    void ParseObject(const std::string& in); 
    void GetObjects(const std::string& in);
    void GetArray  (const std::string& in);
    void Load      (const std::string& in);
    void Print     (int Level=0);
    void Dump      (char* pFile, int Level=0, bool whitespace=true);
    void Dump      (FILE* pFile=stdout, int Level=0);
    std::string DumpToString(int Level=0);
    inline bool   isTag(std::string searchkey) {for(unsigned int i=0;i<key.size();i++){if(key[i] == searchkey)return true;}return false; }
    inline Object& getObject  (std::string searchkey) {for(unsigned int i=0;i<key.size();i++){if(key[i] == searchkey){return obj[i];}}  key.push_back(searchkey); obj.push_back(Object()); return obj[obj.size()-1]; }
    inline Object& operator[] (std::string searchkey) {return getObject(searchkey); }
    inline Object& operator[] (int i                ) {return obj[i]; }
    inline const char* c_str(){return val.c_str();}
    inline std::string toString(){return val;}
    inline double toDouble(){double tmp; sscanf(val.c_str(),"%lf",&tmp);return tmp;}
    inline double toInt   (){int tmp; sscanf(val.c_str(),"%i",&tmp);return tmp;}
    inline bool   toBool  (){if( (val[0]=='t' || val[0]=='T') && (val[1]=='r' || val[1]=='R') && (val[2]=='u' || val[2]=='U') && (val[3]=='e' || val[3]=='E') )return true; return false; }
    inline std::string getString(std::string searchkey, std::string defaultValue=""   ){if(!isTag(searchkey)){return defaultValue;}else{ return getObject(searchkey).toString();  } }
    inline std::string getFullString(std::string searchkey, std::string defaultValue=""   ){if(!isTag(searchkey)){return defaultValue;}else{ return getObject(searchkey).DumpToString();  } }
    inline double getDouble(std::string searchkey, double      defaultValue=0.0  ){if(!isTag(searchkey)){return defaultValue;}else{ return getObject(searchkey).toDouble();  } }
    inline int    getInt   (std::string searchkey, int         defaultValue=0    ){if(!isTag(searchkey)){return defaultValue;}else{ return getObject(searchkey).toInt   ();  } }
    inline bool   getBool  (std::string searchkey, bool        defaultValue=false){if(!isTag(searchkey)){return defaultValue;}else{ return getObject(searchkey).toBool  ();  } }
    inline std::vector<Object>& daughters(){return obj;}
    inline void add(std::string newkey, std::string newval, int length=-1){key.push_back(newkey); obj.push_back(Object(newval, false, length));}
    inline void add(std::string newkey, double newval){char buffer[255];sprintf(buffer,"%f",newval); add(newkey,buffer);}
    inline void addList (){key.push_back("obj"); obj.push_back(Object("LIST"));}
    inline void addArray(std::string name){key.push_back(name); obj.push_back(Object("ARRAY"));}
    inline void setValue(std::string val_){val=val_;}

    inline bool isNumber(){for(std::string::iterator it=val.begin(); it!=val.end();it++){if((*it)<43 || (*it)>57)return false;} return true;}
    inline bool isBool  (){if(val=="true" || val=="TRUE" || val=="True" || val=="false" || val=="FALSE" || val=="False")return true; return false;}
    inline bool isString(){return !isNumber() && !isBool();}

    int EndOfObject;
    std::vector<std::string> key;
    std::vector<Object> obj;
    std::string val;
    int valLength;
  };

}


#endif
