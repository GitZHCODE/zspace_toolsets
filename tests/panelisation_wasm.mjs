import assert from 'node:assert/strict';
import fs from 'node:fs';
const factory=(await import('../wasm/out/zspace_toolsets.js')).default;
const module=await factory();
const call=(name,args=[],types=args.map(()=>'number'))=>module.ccall(name,'number',types,args);
const error=()=>module.UTF8ToString(call('zspace_toolsets_last_error'));
const handle=call('zspace_panelisation_create');assert.ok(handle);
const compute=input=>{module.FS.writeFile('/panels.json',JSON.stringify(input));assert.equal(call('zspace_panelisation_compute',[handle,'/panels.json'],['number','string']),1,error());return JSON.parse(module.UTF8ToString(call('zspace_panelisation_result',[handle])));};
const ribbon={vertexPositions:[0,0,0,0,1,0,0,2,0,2,0,0,2,1,0,2,2,0,2,0,2,2,1,2,2,2,2],faceCounts:[4,4,4,4],faceConnects:[0,3,4,1,1,4,5,2,3,6,7,4,4,7,8,5],creaseEdges:[1,4,4,7]};
const input={meshes:[ribbon],settings:{spacing:1,axis:1}};
const mesh=()=>{
  assert.equal(call('zspace_panelisation_select_mesh',[handle]),1,error());
  const copy=(name,count,float=false)=>{const offset=call(name,[handle])/4,length=call(count,[handle]);return Array.from((float?module.HEAPF32:module.HEAPU32).slice(offset,offset+length));};
  return {vertexPositions:copy('zspace_toolsets_positions','zspace_toolsets_positions_count',true),faceCounts:copy('zspace_toolsets_face_counts','zspace_toolsets_face_count'),faceConnects:copy('zspace_toolsets_face_connects','zspace_toolsets_face_connects_count')};
};
const area=m=>{let cursor=0,sum=0;const p=m.vertexPositions;for(const n of m.faceCounts){const f=m.faceConnects.slice(cursor,cursor+n);cursor+=n;for(let j=1;j+1<n;j++){const a=f[0]*3,b=f[j]*3,c=f[j+1]*3;const u=[0,1,2].map(k=>p[b+k]-p[a+k]),v=[0,1,2].map(k=>p[c+k]-p[a+k]);sum+=Math.hypot(u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])/2;}}return sum;};
let result=compute(input);assert.equal(result.panels.length,8);assert.ok(result.panels.every(p=>Math.abs(p.area-1)<1e-6));assert.ok(Math.abs(area(mesh())-8)<1e-6);
for(const reverse of [false,true]){result=compute({...input,settings:{spacing:0.75,offset:0.2,axis:1,reverse}});assert.ok(Math.abs(area(mesh())-8)<1e-5);assert.ok(result.jointKinds.includes(1));assert.ok(result.jointKinds.includes(2));}
// Separate Rhino objects share one seam; cuts must reuse its vertices.
const left={vertexPositions:[0,0,0,2,0,0,2,1,0,0,1,0],faceCounts:[4],faceConnects:[0,1,2,3],creaseEdges:[]};
const right={vertexPositions:[0,1,0,2,1,0,2,2,0,0,2,0],faceCounts:[4],faceConnects:[0,1,2,3],creaseEdges:[]};
result=compute({meshes:[left,right],settings:{axis:1,spacing:0.5}});assert.equal(result.weldedEdges,1);assert.equal(result.components,1);assert.equal(result.panels.length,8);
const seamMesh=mesh(),seamVertices=[];for(let i=0;i<seamMesh.vertexPositions.length;i+=3)if(Math.abs(seamMesh.vertexPositions[i+1]-1)<1e-6)seamVertices.push(seamMesh.vertexPositions[i]);assert.equal(seamVertices.filter(x=>Math.abs(x-0.5)<1e-6).length,1);
// Independent dimensions on a smooth arch with no world-axis stationing.
const arch={vertexPositions:[],faceCounts:[],faceConnects:[],creaseEdges:[]};
for(let i=0;i<=16;i++)for(let j=0;j<3;j++)arch.vertexPositions.push(2*Math.cos(i*Math.PI/16),j,2*Math.sin(i*Math.PI/16));
for(let i=0;i<16;i++)for(let j=0;j<2;j++){const a=i*3+j;arch.faceCounts.push(4);arch.faceConnects.push(a,a+3,a+4,a+1);}
result=compute({meshes:[arch],settings:{spacing:0.5,width:0.5}});
assert.equal(result.panels.length,52);assert.ok(result.jointKinds.includes(3));assert.ok(Math.abs(area(mesh())-area(arch))<1e-5);
// A mixed triangle/quad wall must form a regular grid, without diagonal wavefronts.
const mixedWall={vertexPositions:[0,0,0,1,0,0,3,0,0,0,2,0,1,2,0,3,2,0],faceCounts:[4,3,3],faceConnects:[0,1,4,3,1,2,5,1,5,4],creaseEdges:[]};
result=compute({meshes:[mixedWall],settings:{spacing:1,width:0.5}});
assert.equal(result.panels.length,12);assert.ok(result.panels.every(p=>Math.abs(p.area-0.5)<1e-5));
const wallMesh=mesh();for(let i=0;i<wallMesh.vertexPositions.length;i+=3){const x=wallMesh.vertexPositions[i],y=wallMesh.vertexPositions[i+1];assert.ok(Math.abs(x-Math.round(x))<1e-5||Math.abs(y*2-Math.round(y*2))<1e-5,'Cuts align with intrinsic wall rows');}
// Patch-specific wall frame: centred vertical cuts and common horizontal heights.
const upright={vertexPositions:[-2,0,0,2,0,0,2,0,3,-2,0,3],faceCounts:[4],faceConnects:[0,1,2,3],creaseEdges:[],panelFrame:[1,0,0,0,0,0,1,0]};
result=compute({meshes:[upright],settings:{spacing:1,width:0.75}});
assert.equal(result.panels.length,16);assert.ok(result.panels.every(p=>Math.abs(p.area-0.75)<1e-6));
assert.deepEqual([...new Set(result.panels.map(p=>p.station))].sort((a,b)=>a-b),[-2,-1,0,1]);
// A crease ending inside a chart anchors the top; do not leave a 0.05 m sliver
// between its x=0.95 rail and an independent x=1 equal-spacing cut.
const creasedWall={vertexPositions:[],faceCounts:[],faceConnects:[],creaseEdges:[1,5],topologyGuides:true};
for(const z of [2,1,0])for(const x of [0,0.95,2,3])creasedWall.vertexPositions.push(x,0,z);
for(let row=0;row<2;row++)for(let col=0;col<3;col++){const a=row*4+col;creasedWall.faceCounts.push(4);creasedWall.faceConnects.push(a,a+1,a+5,a+4);}
const topStarts=r=>[...new Set(r.guides.positions.reduce((out,_,i,values)=>{if(i%3===0&&Math.abs(values[i+2]-2)<1e-6)out.push(Number(values[i].toFixed(5)));return out;},[]))].sort((a,b)=>a-b);
result=compute({meshes:[creasedWall],settings:{spacing:1,width:1}});
const anchored=topStarts(result);assert.equal(anchored.length,5);
for(let i=0;i<5;i++)assert.ok(Math.abs(anchored[i]-[0,.95,.95+2.05/3,.95+4.1/3,3][i])<1e-5,'Equal spacing is redistributed between crease anchors');
assert.ok(Math.abs(area(mesh())-6)<1e-5);
const creaseOutput=mesh();let continuedRedLength=0;
for(let i=0;i<result.jointKinds.length;i++){
 const [a,b]=result.jointEdges.slice(i*2,i*2+2).map(v=>creaseOutput.vertexPositions.slice(v*3,v*3+3));
 if(Math.abs(a[0]-.95)<1e-5&&Math.abs(b[0]-.95)<1e-5&&Math.max(a[2],b[2])<=1+1e-5){
  assert.equal(result.jointKinds[i],2,'Crease continuation stays red after the authored crease ends');continuedRedLength+=Math.abs(a[2]-b[2]);
 }
}
assert.ok(Math.abs(continuedRedLength-1)<1e-5);

result=compute({meshes:[{...creasedWall,creaseEdges:[]}],settings:{spacing:1,width:1}});
assert.deepEqual(topStarts(result),[0,1,2,3],'Uncreased top retains equal spacing');

// A widening quad fan must retain its three equal top seeds and insert red paths.
const fan={vertexPositions:[-1,0,2,0,0,2,1,0,2,-2,0,1,0,0,1,2,0,1,-3,0,0,0,0,0,3,0,0],faceCounts:[4,4,4,4],faceConnects:[0,1,4,3,1,2,5,4,3,4,7,6,4,5,8,7],creaseEdges:[],topologyGuides:true};
result=compute({meshes:[fan],settings:{spacing:1,width:1}});
assert.equal(result.guides.seeds,3);assert.ok(result.guides.branches>=2);assert.ok(result.guides.maxSampledWidth<=1.000001);assert.equal(result.guides.unsupportedCharts,0);
assert.equal(result.guides.positions.length,result.guides.kinds.length*6);
const fanMesh=mesh();assert.ok(Math.abs(area(fanMesh)-8)<1e-5,'Adaptive panel mesh covers the widening fan');
const edgeIncidence=new Map();let offset=0;for(const n of fanMesh.faceCounts){const f=fanMesh.faceConnects.slice(offset,offset+n);offset+=n;for(let i=0;i<n;i++){const key=[f[i],f[(i+1)%n]].sort((a,b)=>a-b).join(',');edgeIncidence.set(key,(edgeIncidence.get(key)||0)+1);}}
for(const [key,count] of edgeIncidence){assert.ok(count<=2);if(count===1){const [a,b]=key.split(',').map(Number).map(i=>fanMesh.vertexPositions.slice(i*3,i*3+3));assert.ok([0,2].some(z=>Math.abs(a[2]-z)<1e-5&&Math.abs(b[2]-z)<1e-5)||[-1,1].some(sign=>Math.abs(a[0]-sign*(3-a[2]))<1e-5&&Math.abs(b[0]-sign*(3-b[2]))<1e-5),'No unmatched internal edge at branch/Y junction '+JSON.stringify({a,b}));}}
const coarseFanCount=result.panels.length;
result=compute({meshes:[fan],settings:{spacing:1,width:0.5}});assert.ok(result.panels.length>coarseFanCount,'Smaller Y makes more mesh patches');assert.ok(Math.abs(area(mesh())-8)<1e-5);

// A transverse topology divider must retain the flow field and stop Y joints at the preceding flow boundary.
const fanWithoutDivider = compute({meshes:[fan],settings:{spacing:1,width:0.4}});
result=compute({meshes:[{...fan,panelSeamSeeds:[3,4]}],settings:{spacing:1,width:0.4}});
assert.deepEqual(result.guides,fanWithoutDivider.guides,'Divider does not restart flow stationing');
assert.ok(result.jointKinds.includes(3),'Y joints remain away from the divider');
assert.ok(Math.abs(area(mesh())-8)<1e-5);
function checkDivider(result,m,expectedLength) {
  const incidence=new Map();let cursor=0;
  for(let fi=0;fi<m.faceCounts.length;fi++) {const n=m.faceCounts[fi],f=m.faceConnects.slice(cursor,cursor+n);cursor+=n;
    for(let j=0;j<n;j++){const key=[f[j],f[(j+1)%n]].sort((a,b)=>a-b).join(',');if(!incidence.has(key))incidence.set(key,[]);incidence.get(key).push(fi);}
  }
  let length=0;const seamVertices=new Set(),yVertices=new Set(),degree=new Map();
  for(let j=0;j<result.jointKinds.length;j++){
    const [a,b]=result.jointEdges.slice(j*2,j*2+2),kind=result.jointKinds[j];
    if(kind===3){yVertices.add(a);yVertices.add(b);}
    if(kind!==4)continue;
    const faces=incidence.get([a,b].sort((a,b)=>a-b).join(','));
    assert.equal(faces.length,2,'Every divider segment belongs to both neighbouring mesh patches');
    assert.notEqual(result.facePanels[faces[0]],result.facePanels[faces[1]]);
    for(const v of [a,b]){seamVertices.add(v);degree.set(v,(degree.get(v)||0)+1);}
    length+=Math.hypot(...[0,1,2].map(k=>m.vertexPositions[a*3+k]-m.vertexPositions[b*3+k]));
  }
  assert.ok(Math.abs(length-expectedLength)<1e-4,'Complete curved rail preserved');
  assert.equal([...degree.values()].filter(n=>n===1).length,2,'Continuous divider has only two ends');
  assert.ok([...degree.values()].every(n=>n<=2),'Divider does not branch');
  assert.equal([...seamVertices].filter(v=>yVertices.has(v)).length,0,'Y joints stop before the curved divider');
}
checkDivider(result,mesh(),4);

// Arch rails share an odd division count; adjacent strips use alternate indices.
const loopSheet={vertexPositions:[],faceCounts:[],faceConnects:[],creaseEdges:[],layoutMode:2,staggerFraction:0.5};
for(let y=0;y<=2;y++)for(let x=0;x<=16;x++)loopSheet.vertexPositions.push(x/4,y,0);
for(let y=0;y<2;y++)for(let x=0;x<16;x++){const a=y*17+x;loopSheet.faceCounts.push(4);loopSheet.faceConnects.push(a,a+1,a+18,a+17);}
for(let x=0;x<16;x++)loopSheet.creaseEdges.push(17+x,18+x);
result=compute({meshes:[loopSheet],settings:{spacing:1,width:2}});const loopMesh=mesh();
assert.equal(result.panels.length,4);assert.ok(Math.abs(area(loopMesh)-8)<1e-5);
const loopEdges=new Map();let lc=0;for(const n of loopMesh.faceCounts){const f=loopMesh.faceConnects.slice(lc,lc+n);lc+=n;for(let j=0;j<n;j++){const key=[f[j],f[(j+1)%n]].sort((a,b)=>a-b).join(',');loopEdges.set(key,(loopEdges.get(key)||0)+1);}}
for(const [key,count] of loopEdges){assert.ok(count<=2);if(count===1){const points=key.split(',').map(Number).map(v=>loopMesh.vertexPositions.slice(v*3,v*3+3));assert.ok([[0,0],[0,4],[1,0],[1,2]].some(([axis,value])=>points.every(p=>Math.abs(p[axis]-value)<1e-6)),'Staggered neighbours share seam vertices without internal cracks');}}

const endCuts=y=>{const xs=[];for(let i=0;i<result.jointKinds.length;i++)if(result.jointKinds[i]===1)for(const v of result.jointEdges.slice(i*2,i*2+2)){const p=loopMesh.vertexPositions.slice(v*3,v*3+3);if(Math.abs(p[1]-y)<1e-6)xs.push(p[0]);}return [...new Set(xs)].sort((a,b)=>a-b);};
assert.equal(endCuts(0).length,1);assert.ok(Math.abs(endCuts(0)[0]-8/3)<1e-5);
assert.equal(endCuts(2).length,1);assert.ok(Math.abs(endCuts(2)[0]-4/3)<1e-5);
// ceil(4/(2*.5))=4 is raised to the next odd count, 5.
const odd=compute({meshes:[loopSheet],settings:{spacing:.5,width:2}}), oddMesh=mesh();
assert.equal(odd.panels.length,6);assert.ok(Math.abs(area(oddMesh)-8)<1e-5);
for(let i=0;i<odd.jointKinds.length;i++)if(odd.jointKinds[i]===1)for(const v of odd.jointEdges.slice(i*2,i*2+2))assert.ok(Math.abs(oddMesh.vertexPositions[v*3]/.8-Math.round(oddMesh.vertexPositions[v*3]/.8))<1e-5,'Arch cuts use five equal rail divisions');

// Six-metre rails: n=ceil(6/(2*1))=3; alternate rows use even and odd divisions at 4m / 2m.
const alignWall={...creasedWall,vertexPositions:creasedWall.vertexPositions.map((v,i)=>i%3===0?v+1:v)};
const ceiling={vertexPositions:[],faceCounts:[],faceConnects:[],creaseEdges:[],layoutMode:1,alignmentSource:0,alignmentDirection:[1,0,0]};
for(let y=0;y<=2;y++)for(let x=0;x<=24;x++)ceiling.vertexPositions.push(x/4,y,2);
for(let y=0;y<2;y++)for(let x=0;x<24;x++){const a=y*25+x;ceiling.faceCounts.push(4);ceiling.faceConnects.push(a,a+1,a+26,a+25);}
result=compute({meshes:[alignWall,ceiling],settings:{spacing:1,width:1}});const ceilingMesh=mesh();
assert.ok(Math.abs(area(ceilingMesh)-18)<1e-5);
const ceilingCutXs=[];for(let i=0;i<result.jointKinds.length;i++)if(result.jointKinds[i]===1){const p=result.jointEdges.slice(i*2,i*2+2).map(v=>ceilingMesh.vertexPositions.slice(v*3,v*3+3));if(p.every(q=>Math.abs(q[2]-2)<1e-6)&&Math.abs(p[0][0]-p[1][0])<1e-6)ceilingCutXs.push(p[0][0]);}
for(const x of [2,4])assert.ok(ceilingCutXs.some(v=>Math.abs(v-x)<1e-5),'Ceiling alternate equal rail division '+x);
for(let i=0;i<result.jointKinds.length;i++)if(result.jointKinds[i]===1){const ps=result.jointEdges.slice(i*2,i*2+2).map(v=>ceilingMesh.vertexPositions.slice(v*3,v*3+3));
  if(ps.every(p=>Math.abs(p[2]-2)<1e-6))for(const p of ps){if(Math.abs(p[1])<1e-6)assert.ok(Math.abs(p[0]-4)<1e-5,'First band uses even divisions');if(Math.abs(p[1]-2)<1e-6)assert.ok(Math.abs(p[0]-2)<1e-5,'Adjacent band uses odd divisions');}}


// Offset rail endpoints must not incline cross-cuts on a straight reach.
const skewedCeiling={...ceiling,alignmentSource:0,vertexPositions:ceiling.vertexPositions.map((v,i)=>i%3===0?v+0.2*ceiling.vertexPositions[i+1]:v)};
result=compute({meshes:[skewedCeiling],settings:{spacing:1,width:1}});const skewedMesh=mesh();let straightCuts=0;
for(let j=0;j<result.jointKinds.length;j++)if(result.jointKinds[j]===1){const [a,b]=result.jointEdges.slice(j*2,j*2+2).map(v=>skewedMesh.vertexPositions.slice(v*3,v*3+3));
  assert.ok(Math.abs(a[0]-b[0])<1e-6,'Straight ceiling cross-cuts remain parallel despite offset rail starts');straightCuts++;}
assert.ok(straightCuts>0);assert.ok(Math.abs(area(skewedMesh)-12)<1e-5);

// Four rows use alternating A,B,A,B stations.
const pairedRoof={vertexPositions:[],faceCounts:[],faceConnects:[],layoutMode:1,alignmentSource:0,alignmentDirection:[1,0,0]};
for(let y=0;y<=4;y++)for(let x=0;x<=24;x++)pairedRoof.vertexPositions.push(x/4,y,2);
for(let y=0;y<4;y++)for(let x=0;x<24;x++){const a=y*25+x;pairedRoof.faceCounts.push(4);pairedRoof.faceConnects.push(a,a+1,a+26,a+25);}
result=compute({meshes:[pairedRoof],settings:{spacing:1,width:1}});const pairedMesh=mesh();
for(let j=0;j<result.jointKinds.length;j++)if(result.jointKinds[j]===1){const ps=result.jointEdges.slice(j*2,j*2+2).map(v=>pairedMesh.vertexPositions.slice(v*3,v*3+3)),y=(ps[0][1]+ps[1][1])/2;
  if(Math.abs(ps[0][1]-ps[1][1])>1e-6)for(const p of ps)assert.ok(Math.abs(p[0]-(Math.floor(y)%2===0?4:2))<1e-5,'A,B,A,B ceiling row staggering');}
assert.ok(Math.abs(area(pairedMesh)-24)<1e-5);

const before=JSON.stringify(result);module.FS.writeFile('/panels.json',JSON.stringify({...input,settings:{spacing:0}}));assert.equal(call('zspace_panelisation_compute',[handle,'/panels.json'],['number','string']),0);assert.match(error(),/spacing/);assert.deepEqual(JSON.parse(module.UTF8ToString(call('zspace_panelisation_result',[handle]))),JSON.parse(before));
const dataPath=new URL('../../zspace_alice_webviewer/sketches/panelisationInput.json',import.meta.url);
if(fs.existsSync(dataPath)){
  const {preparePanelisationInput}=await import('../../zspace_alice_webviewer/sketches/panelisationInput.js');
  const data=JSON.parse(fs.readFileSync(dataPath,'utf8'));const {sources,...prepared}=preparePanelisationInput(data,{});
  const baseline=compute(preparePanelisationInput(data,{roofReview:false}));
  result=compute(prepared);const output=mesh();
  for(const id of ['d91c1b65-62c7-45fc-9da9-639782c4c251','387cbe7b-dda6-4966-a667-78c17bdd7710','95c3b292-a6f0-4116-aac6-4fafe8c4b1f0']) {
    const source=sources.findIndex(o=>o.id===id);
    assert.deepEqual(result.panels.filter(p=>p.source===source).map(({chart,...p})=>p),baseline.panels.filter(p=>p.source===source).map(({chart,...p})=>p),'Approved wall panel layout unchanged');
  }
  // Each custom transverse joint lies in a single section plane, including
  // alternate staggered strips. Test the actual panel mesh edges, not overlays.
  const edgeFaces=new Map();let offset=0;
  for(let f=0;f<output.faceCounts.length;f++){const n=output.faceCounts[f],ids=output.faceConnects.slice(offset,offset+n);offset+=n;
    for(let j=0;j<n;j++){const key=[ids[j],ids[(j+1)%n]].sort((a,b)=>a-b).join(',');if(!edgeFaces.has(key))edgeFaces.set(key,[]);edgeFaces.get(key).push(f);}}
  const cuts=new Map();
  for(let j=0;j<result.jointKinds.length;j++)if(result.jointKinds[j]===1){const ids=result.jointEdges.slice(2*j,2*j+2),faces=edgeFaces.get([...ids].sort((a,b)=>a-b).join(','))||[];
    if(faces.length!==2)continue;const source=result.faceSources[faces[0]],input=prepared.meshes[source];if(!input.layoutMode)continue;
    const key=source+':'+faces.map(f=>result.facePanels[f]).sort((a,b)=>a-b).join(',');if(!cuts.has(key))cuts.set(key,{source,ids:new Set()});for(const v of ids)cuts.get(key).ids.add(v);}
  assert.ok(cuts.size>100,'Custom layout has measurable transverse panel joints');
  for(const {source,ids} of cuts.values()){const input=prepared.meshes[source],h=input.alignmentDirection;
    if(input.layoutMode===1){const points=[...ids].map(v=>output.vertexPositions.slice(3*v,3*v+3)),a=points[0];
      const sub=(p,q)=>p.map((v,k)=>v-q[k]),cross=(u,v)=>[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];
      const b=points.reduce((best,p)=>Math.hypot(...sub(p,a))>Math.hypot(...sub(best,a))?p:best,a),u=sub(b,a);
      const normal=points.map(p=>cross(u,sub(p,a))).reduce((best,n)=>Math.hypot(...n)>Math.hypot(...best)?n:best,[0,0,0]),length=Math.hypot(...normal);
      if(length>1e-8)for(const p of points)assert.ok(Math.abs(sub(p,a).reduce((sum,v,k)=>sum+v*normal[k],0))/length<2e-5,'Ceiling cuts stay in a section plane through the ground return');
      continue;
    }
    const points=[...ids].map(v=>output.vertexPositions.slice(3*v,3*v+3)).map(p=>input.layoutMode===1?[p[0],p[1]]:[p[0]*h[0]+p[1]*h[1],p[2]]);
    const a=points[0],b=points.reduce((best,p)=>Math.hypot(p[0]-a[0],p[1]-a[1])>Math.hypot(best[0]-a[0],best[1]-a[1])?p:best,a),dx=b[0]-a[0],dy=b[1]-a[1],length=Math.hypot(dx,dy);
    if(length>1e-5)for(const p of points)assert.ok(Math.abs(dx*(p[1]-a[1])-dy*(p[0]-a[0]))/length<2e-5,'Clean section cut on source '+source);
  }
  checkDivider(result,output,19.37751699);
  assert.equal(prepared.meshes.filter(m=>m.panelSeamSeeds).length,1);
  assert.equal(prepared.meshes.filter(m=>m.topologyGuides).length,3,'Only the three reviewed wall/fillet patches use the new method');assert.ok(output.vertexPositions.every(Number.isFinite));assert.equal(new Set(result.faceSources).size,sources.length);
  assert.ok(result.guides.kinds.includes(0));assert.ok(result.guides.kinds.includes(1));assert.ok(result.guides.maxSampledWidth<=prepared.settings.spacing*1.00001);

  console.log({guideSeeds:result.guides.seeds,branches:result.guides.branches,unsupported:result.guides.unsupportedCharts});
  const inputArea=prepared.meshes.reduce((n,m)=>n+area(m),0),outputArea=area(output);
  assert.ok(Math.abs(inputArea-outputArea)<inputArea*1e-5,'Area conserved across all source walls');
  assert.equal(result.facePanels.length,output.faceCounts.length);assert.ok(result.facePanels.every(i=>i>=0&&i<result.panels.length));
  console.log(JSON.stringify({walls:sources.length,panels:result.panels.length,components:result.components,weldedEdges:result.weldedEdges,smallPanels:result.smallPanels,inputArea,outputArea}));
  const review = compute(preparePanelisationInput(data,{spacing:1,width:2})), reviewMesh=mesh();
  checkDivider(review,reviewMesh,19.37751699);
  const ceilingSource=sources.findIndex(o=>o.id==='0e9e6737-d96e-4ed5-8e26-7821be8a036d');
  const ceilingChartCounts=new Map();for(const p of review.panels)if(p.source===ceilingSource)ceilingChartCounts.set(p.chart,(ceilingChartCounts.get(p.chart)||0)+1);
  const counts=[...ceilingChartCounts.values()].sort((a,b)=>b-a);
  assert.ok(counts.length>1,'Fixture includes extraordinary end patches');
  assert.ok(counts.slice(1).every(n=>n<=8),'Short end patches inherit main stations instead of restarting all ceiling divisions');

  const wallVertices=new Set();let faceCursor=0;
  for(let f=0;f<reviewMesh.faceCounts.length;f++){const n=reviewMesh.faceCounts[f];if(review.faceSources[f]===22)for(const v of reviewMesh.faceConnects.slice(faceCursor,faceCursor+n))wallVertices.add(v);faceCursor+=n;}
  const redVertices=new Set(),cyanDegree=new Map();
  for(let j=0;j<review.jointKinds.length;j++){const vs=review.jointEdges.slice(j*2,j*2+2),kind=review.jointKinds[j];
    if([0,2,4].includes(kind))for(const v of vs)redVertices.add(v);
    if(kind===3&&vs.every(v=>wallVertices.has(v)))for(const v of vs)cyanDegree.set(v,(cyanDegree.get(v)||0)+1);
  }
  assert.ok(cyanDegree.size>0,'Straight wall retains cyan rows');
  for(const [v,degree] of cyanDegree)if(degree===1)assert.ok(redVertices.has(v),'Cyan row ends on a red rail or source boundary, never a blue path');

  if(process.env.ZSPACE_GUIDE_PREVIEW)fs.writeFileSync(process.env.ZSPACE_GUIDE_PREVIEW,JSON.stringify({mesh:reviewMesh,...review}));
}
assert.equal(call('zspace_toolsets_destroy',[handle]),1);assert.equal(call('zspace_panelisation_select_mesh',[handle]),0);
console.log('Panelisation WASM: wall rails, ceiling alignment, topology-snapped stagger, shared seams and area conservation passed.');
