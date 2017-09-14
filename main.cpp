#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "Gmsh.h"
#include "GModel.h"
#include "MElement.h"
#include "MTriangle.h"
#include "MQuadrangle.h"
#include "MElementCut.h"
#include "MLine.h"
#include "MVertex.h"

enum Graph{
  Node,
  Element
};

struct Param{
  double dx = 0.;
  double dy = 0.;
  double scale = 1.;
  Graph graph = Node;
  bool only = false;
  std::vector<int> onlyList;
};

void writeTikz(GModel *model, Param &param);
template <class ITERATOR>
void fillVtoV(std::unordered_map<MVertex*, std::vector<MVertex*> > &vertexToVertices, Param &param, ITERATOR it_beg, ITERATOR it_end);
template <class ITERATOR>
void fillEtoV(std::unordered_map<MElement*, std::vector<MVertex*> > &elementToVertices, Param &param, ITERATOR it_beg, ITERATOR it_end);
void readParam(char *fileName, Param &param);
void readList(std::vector<int> &onlyList, std::string &second);

int main(int argc, char **argv)
{
  GmshInitialize(1, argv);
  GModel *model = new GModel();
  std::cout << "Reading msh file... " << std::flush;
  model->readMSH(argv[1]);
  Param param;
  if(argc > 2)
  {
    readParam(argv[2], param);
  }
  std::cout << "Done!" << std::endl;
  
  unsigned int numElem[6];
  const int meshDim = model->getNumMeshElements(numElem);
  if(meshDim > 2)
  {
    std::cout << "Error: Not implemented for 3D meshes" << std::endl;
    return 1;
  }
  
  writeTikz(model, param);
    
  delete model;
  GmshFinalize();
  return 0;
}

void writeTikz(GModel *model, Param &param)
{
  std::ofstream file("tikz.txt");
  
  if(param.graph == Node)
  {
    std::unordered_map<MVertex*, std::vector<MVertex*> > vertexToVertices;
    
    if(model->getDim() == 2)
    {
      for(GModel::fiter it = model->firstFace(); it != model->lastFace(); ++it)
      {
        fillVtoV(vertexToVertices, param, (*it)->triangles.begin(), (*it)->triangles.end());
        fillVtoV(vertexToVertices, param, (*it)->quadrangles.begin(), (*it)->quadrangles.end());
        fillVtoV(vertexToVertices, param, (*it)->polygons.begin(), (*it)->polygons.end());
      }
    }
    else if(model->getDim() == 1)
    {
      for(GModel::eiter it = model->firstEdge(); it != model->lastEdge(); ++it)
      {
        fillVtoV(vertexToVertices, param, (*it)->lines.begin(), (*it)->lines.end());
      }
    }
    
    for(std::unordered_map<MVertex*, std::vector<MVertex*> >::iterator it = vertexToVertices.begin(); it != vertexToVertices.end(); ++it)
    {
      MVertex* v0 = (*it).first;
      std::vector<MVertex*> vertices = (*it).second;
      
      for(unsigned int i = 0; i < vertices.size(); i++)
      {
        MVertex* v1 = vertices[i];
        
        file << "\\draw (" << param.scale*v0->x()+param.dx << "," << param.scale*v0->y()+param.dy << ") -- (" << param.scale*v1->x()+param.dx << "," << param.scale*v1->y()+param.dy << ");" << std::endl;
      }
    }
  }
  else if(param.graph == Element)
  {
    std::unordered_map<MElement*, std::vector<MVertex*> > elementToVertices;
    
    if(model->getDim() == 2)
    {
      for(GModel::fiter it = model->firstFace(); it != model->lastFace(); ++it)
      {
        fillEtoV(elementToVertices, param, (*it)->triangles.begin(), (*it)->triangles.end());
        fillEtoV(elementToVertices, param, (*it)->quadrangles.begin(), (*it)->quadrangles.end());
        fillEtoV(elementToVertices, param, (*it)->polygons.begin(), (*it)->polygons.end());
      }
    }
    else if(model->getDim() == 1)
    {
      for(GModel::eiter it = model->firstEdge(); it != model->lastEdge(); ++it)
      {
        fillEtoV(elementToVertices, param, (*it)->lines.begin(), (*it)->lines.end());
      }
    }
    
    for(std::unordered_map<MElement*, std::vector<MVertex*> >::iterator it = elementToVertices.begin(); it != elementToVertices.end(); ++it)
    {
      std::vector<MVertex*> vertices = (*it).second;
      file << "\\draw[fill = white] ";
      for(unsigned int i = 0; i < vertices.size(); i++)
      {
         file << "(" << param.scale*vertices[i]->x()+param.dx << "," << param.scale*vertices[i]->y()+param.dy << ") -- ";
      }
      file << " cycle;" << std::endl;
    }
  }
  
  file.close();
}

template <class ITERATOR>
void fillVtoV(std::unordered_map<MVertex*, std::vector<MVertex*> > &vertexToVertices, Param &param, ITERATOR it_beg, ITERATOR it_end)
{
  for(ITERATOR IT = it_beg; IT != it_end ; ++IT)
  {
    MElement *el = *IT;
    
    if(param.only)
    {
      std::vector<int>::iterator it = find(param.onlyList.begin(), param.onlyList.end(), el->getNum());
      if(it == param.onlyList.end()) continue;
    }
    
    for(unsigned int j = 0; j < el->getNumEdges(); j++)
    {
      const MEdge e = el->getEdge(j);
      if(e.getVertex(0)->getNum() <= e.getVertex(1)->getNum())
      {
        bool found = false;
        for(unsigned int k = 0; k < vertexToVertices[e.getVertex(0)].size(); k++)
        {
          if(vertexToVertices[e.getVertex(0)][k] == e.getVertex(1))
          {
            found = true;
            break;
          }
        }
                
        if(!found) vertexToVertices[e.getVertex(0)].push_back(e.getVertex(1));
      }
      else
      {
        bool found = false;
        for(unsigned int k = 0; k < vertexToVertices[e.getVertex(1)].size(); k++)
        {
          if(vertexToVertices[e.getVertex(1)][k] == e.getVertex(0))
          {
            found = true;
            break;
          }
        }
        
        if(!found) vertexToVertices[e.getVertex(1)].push_back(e.getVertex(0));
      }
    }
  }
}

template <class ITERATOR>
void fillEtoV(std::unordered_map<MElement*, std::vector<MVertex*> > &elementToVertices, Param &param, ITERATOR it_beg, ITERATOR it_end)
{
  for(ITERATOR IT = it_beg; IT != it_end ; ++IT)
  {
    MElement *el = *IT;
    
    if(param.only)
    {
      std::vector<int>::iterator it = find(param.onlyList.begin(), param.onlyList.end(), el->getNum());
      if(it != param.onlyList.end()) continue;
    }
    
    std::vector<MVertex*> vertices;
    for(unsigned int j = 0; j < el->getNumVertices(); j++)
    {
      vertices.push_back(el->getVertex(j));
    }
    elementToVertices.insert( std::pair<MElement*, std::vector<MVertex*> >(el,vertices) );
  }
}

void readParam(char *fileName, Param &param)
{
  std::ifstream file(fileName);
  std::string line;
  std::string first;
  std::string second;
  
  while(!file.eof())
  {
    file >> line;
    std::size_t pos = line.find("=");
    first = line.substr(0,pos);
    second = line.substr(pos+1);
    
    if(first == "dx")
    {
      param.dx = std::stod(second);
    }
    else if(first == "dy")
    {
      param.dy = std::stod(second);
    }
    else if(first == "graph")
    {
      if(second == "node")
      {
        param.graph = Node;
      }
      else if(second == "element")
      {
        param.graph = Element;
      }
    }
    else if(first == "scale")
    {
      param.scale = std::stod(second);
    }
    else if(first == "only")
    {
      param.only = true;
      readList(param.onlyList, second);
    }
  }
  
  file.close();
}

void readList(std::vector<int> &onlyList, std::string &second)
{
  int i = 1;
  int nb = 0;
  while(second[i] != '}')
  {
    if(second[i] != ',')
    {
      char c = second[i];
      nb = 10*nb + atoi(&c);
    }
    else
    {
      onlyList.push_back(nb);
      nb = 0;
    }
    
    i++;
  }
  onlyList.push_back(nb);
}
