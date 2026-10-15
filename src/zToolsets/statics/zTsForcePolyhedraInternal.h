#pragma once
#include <zspace/zToolsets/statics/zTsForcePolyhedra.h>
#include <Eigen/Dense>
#include <zspace/zInterface/functionsets/zFnMeshDynamics.h>
#include <zspace/zInterface/functionsets/zFnGraphDynamics.h>

namespace zSpace::forcePolyhedra {
using V = Eigen::Vector3d;
inline V vector(const zVector& p) { return V(p.x,p.y,p.z); }
inline zVector publicVector(const V& p) { return zVector(float(p.x()),float(p.y()),float(p.z())); }
inline void require(bool valid, const char* message) { if(!valid) throw std::invalid_argument(message); }
template<class Object,class FunctionSet> struct DynamicsState {
    Object object;
    FunctionSet fn;
    double scale=1;
    DynamicsState():fn(object){}
};
using MeshDynamicsState=DynamicsState<zObjectMesh,zFnMeshDynamics>;
using GraphDynamicsState=DynamicsState<zObjectGraph,zFnGraphDynamics>;
// Integrate small increments in rebased local coordinates. This retains double
// precision diagram positions while the core particles use float coordinates.
// Overdamped Euler clears velocity/derivatives after each physical time step.
template<class State> std::vector<V> integrateCorrections(State& state,const std::vector<V>& corrections,const zForcePolyhedraSettings& s,double stepLimit){
    auto raw=state.fn.getRawVertexPositions();double largest=0;
    std::vector<V> forces;forces.reserve(corrections.size());
    for(int i=0;i<int(corrections.size());++i){raw[i]=zVector();forces.push_back(corrections[i]*(s.relaxation/state.scale));largest=std::max(largest,forces.back().norm());}
    double limit=stepLimit/(state.scale*s.timeStep*s.timeStep),factor=largest>limit?limit/largest:1;
    for(int i=0;i<int(forces.size());++i){zVector force=publicVector(forces[i]*factor);double length=force.length();if(length>0)state.fn.addLoadForce(length,i,force);}
    state.fn.update(s.timeStep,zEuler,true,true,true);
    std::vector<V> displacement;displacement.reserve(forces.size());raw=state.fn.getRawVertexPositions();
    for(int i=0;i<int(forces.size());++i){V delta=vector(raw[i])*state.scale;require(delta.allFinite(),"Nonfinite dynamics displacement.");displacement.push_back(delta);}
    return displacement;
}
struct Cell {
    zObjectMesh mesh;
    std::vector<V> positions;
    std::vector<zIntArray> polygons;
    std::vector<V> normals;
    std::vector<zForcePolyhedraFace> faces;
    zForcePolyhedraCellReport report;
    std::unique_ptr<MeshDynamicsState> dynamics;
};
void buildCell(Cell& cell, const zForcePolyhedraSettings& settings, bool allowUnbalanced = false);
void pursue(Cell& cell, const zForcePolyhedraSettings& settings, int iterations, bool adaptiveNormals = false);
V faceNormal(const Cell& cell, int face);
void refresh(Cell& cell, const zForcePolyhedraSettings& settings);
}
