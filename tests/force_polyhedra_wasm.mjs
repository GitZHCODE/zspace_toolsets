import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
const factory=(await import('../wasm/out/zspace_toolsets.js')).default;
const m=await factory({wasmBinary:await readFile(new URL('../wasm/out/zspace_toolsets.wasm',import.meta.url))});
const call=(name,args=[],types=args.map(()=>'number'))=>m.ccall(name,'number',types,args);
const error=()=>m.UTF8ToString(call('zspace_forcepolyhedra_last_error'));
const ok=n=>assert.equal(n,1,error());
const result=h=>{const p=call('zspace_forcepolyhedra_result',[h]);assert.ok(p,error());return JSON.parse(m.UTF8ToString(p));};
const compute=(h,input)=>{m.FS.writeFile('/force-polyhedra.json',JSON.stringify(input));return call('zspace_forcepolyhedra_compute',[h,'/force-polyhedra.json'],['number','string']);};
// Independent ZST1 encoder/decoder exercises the same typed contract as Alice.
function encode(input){
  const chunks=[];
  const byte=n=>chunks.push(Uint8Array.of(n));
  const scalar=(n,type='Uint32')=>{const b=new ArrayBuffer(type==='Float64'?8:4);new DataView(b)['set'+type](0,n,true);chunks.push(new Uint8Array(b));};
  const string=s=>{const b=new TextEncoder().encode(s);scalar(b.length);chunks.push(b);};
  function value(v){
    if(v===null){byte(0);return;}if(typeof v==='boolean'){byte(v?2:1);return;}
    if(typeof v==='number'){byte(Number.isSafeInteger(v)?3:4);scalar(v,'Float64');return;}
    if(typeof v==='string'){byte(5);string(v);return;}
    if(Array.isArray(v)){
      const numeric=v.every(x=>typeof x==='number');byte(numeric?9:7);scalar(v.length);
      for(const x of v)numeric?scalar(x,'Float64'):value(x);return;
    }
    byte(6);scalar(Object.keys(v).length);for(const [k,x] of Object.entries(v)){string(k);value(x);}
  }
  scalar(0x3154535a);value(input);const out=new Uint8Array(chunks.reduce((n,b)=>n+b.length,0));let at=0;for(const b of chunks){out.set(b,at);at+=b.length;}return out;
}
function decode(bytes){
  const view=new DataView(bytes.buffer,bytes.byteOffset,bytes.byteLength);let at=0;
  const scalar=(type='Uint32')=>{const n=view['get'+type](at,true);at+=type==='Uint8'?1:type==='Float64'?8:4;return n;};
  const string=()=>{const n=scalar();const s=new TextDecoder().decode(bytes.subarray(at,at+n));at+=n;return s;};
  function value(){
    const tag=scalar('Uint8');if(tag===0)return null;if(tag===1)return false;if(tag===2)return true;
    if(tag===3||tag===4)return scalar('Float64');if(tag===5)return string();
    const n=scalar();if(tag===6){const out={};for(let i=0;i<n;i++){const k=string();out[k]=value();}return out;}
    if(tag===7||tag===8||tag===9)return Array.from({length:n},()=>tag===7?value():scalar(tag===8?'Int32':'Float64'));
    throw Error('Unknown ZST1 tag');
  }
  assert.equal(scalar(),0x3154535a);const out=value();assert.equal(at,bytes.length);return out;
}
function typedCompute(h,input,entry='zspace_forcepolyhedra_compute_from_buffers'){
  const bytes=encode(input),p=call('zspace_alloc',[bytes.length]);assert.ok(p);
  try{new Uint8Array(m.HEAPU32.buffer,p,bytes.length).set(bytes);return call(entry,[h,p,bytes.length]);}
  finally{call('zspace_free',[p]);}
}
function typedResult(h){
  const p=call('zspace_forcepolyhedra_result_buffers',[h]);assert.ok(p,error());
  const n=call('zspace_forcepolyhedra_result_buffer_size',[h]);return decode(new Uint8Array(m.HEAPU32.buffer,p,n).slice());
}
const h=call('zspace_forcepolyhedra_create'),other=call('zspace_forcepolyhedra_create');
const arch={graph:{vertexPositions:[-1,0,0,0,0,1,1,0,0],edgeConnects:[0,1,1,2]},supports:[0,2],masses:[1,2,1],edgeForceMagnitudes:[Math.SQRT2,Math.SQRT2]};
try{
  ok(compute(h,arch));const r=result(h);assert.ok(r.equilibrium&&r.cellsConverged);
  assert.ok(Math.abs(r.axialForces[0]+Math.SQRT2)<1e-9);assert.ok(Math.abs(r.supportReactions[2]-2)<1e-6);
  assert.equal(r.cells[1].forceRank,2);assert.equal(r.cells[1].faces.filter(f=>f.type===4).length,2);
  for(const cell of r.cells){assert.ok(cell.vertexPositions.every(Number.isFinite));assert.equal(cell.faceCounts.length,cell.faces.length);}
  ok(typedCompute(other,arch));assert.deepEqual(typedResult(other),r);assert.deepEqual(result(other),r);
  const members=[];
  for(const cell of r.cells)for(const face of cell.faces)for(const source of face.sources)if(source.type===0)members.push({node:cell.vertexId,id:source.id,magnitude:source.magnitude,normal:face.normal});
  for(let e=0;e<2;e++){
    const pair=members.filter(f=>f.id===e);assert.equal(pair.length,2);assert.equal(pair[0].magnitude,pair[1].magnitude);
    for(let k=0;k<3;k++)assert.ok(Math.abs(pair[0].normal[k]+pair[1].normal[k])<1e-6);
  }
  for(const bad of [{...arch,masses:[-1,2,1]},{...arch,supports:[0,-1]},
    {...arch,graph:{...arch.graph,edgeConnects:[0,99]}},{...arch,forceDensities:[0,0]},
    {...arch,settings:{forceScale:0}},{...arch,settings:{timeStep:0}},{...arch,settings:{maxVertexStepRatio:1}},{...arch,settings:{normalStepDegrees:0}},{...arch,loads:[0,0,0]},
    {...arch,graph:{...arch.graph,edgeConnects:[0,1.5]}},{...arch,supports:[4294967296]}]){
    assert.equal(typedCompute(other,bad),0);assert.ok(error());assert.deepEqual(result(other),r);
  }
  assert.equal(call('zspace_forcepolyhedra_compute_from_buffers',[other,0,16]),0);assert.match(error(),/memory/);
  ok(typedCompute(other,{vertexPositions:[-1,0,0,.2,0,.7,1,0,0]},'zspace_forcepolyhedra_set_form_from_buffers'));
  ok(call('zspace_forcepolyhedra_set_form_weight',[other,1]));const forceBefore=result(other).cells;
  ok(call('zspace_forcepolyhedra_equilibrium',[other,2000]));assert.ok(result(other).equilibrium);assert.deepEqual(result(other).cells,forceBefore);
  assert.deepEqual(result(h),r);
  ok(compute(other,{...arch,supports:[]}));assert.equal(result(other).equilibrium,false);assert.equal(result(other).cells.length,3);
  assert.equal(call('zspace_forcepolyhedra_equilibrium',[other,0]),0);
  const branch={graph:{vertexPositions:[1,0,0,-.5,.8660254,0,-.5,-.8660254,0,0,0,1],edgeConnects:[0,3,1,3,2,3]},supports:[0,1,2],masses:[0,0,0,3],edgeForceMagnitudes:[Math.SQRT2,Math.SQRT2,Math.SQRT2]};
  ok(typedCompute(other,branch));let b=typedResult(other);assert.ok(b.equilibrium&&b.cellsConverged);assert.equal(b.cells[3].forceRank,3);
  ok(typedCompute(other,{...branch,settings:{maxIterations:1}}));b=typedResult(other);assert.ok(!b.cellsConverged&&!b.equilibrium);
  ok(call('zspace_forcepolyhedra_step',[other,5000]));assert.ok(typedResult(other).cellsConverged);
  const skew={graph:{...arch.graph,vertexPositions:[-1,0,0,.3,0,.6,1,0,0]},supports:[0,2],masses:[0,2,0],settings:{formWeight:.5},solveEquilibrium:true};
  ok(typedCompute(other,skew));b=typedResult(other);assert.ok(b.equilibrium);assert.ok(b.maxAngleDeviation<.01);assert.ok(b.maxResidual<1e-6);
  assert.equal(call('zspace_forcepolyhedra_set_form_weight',[other,1.1]),0);assert.deepEqual(typedResult(other),b);
  const fixed={...skew,settings:{formWeight:0}};
  ok(typedCompute(other,fixed));b=typedResult(other);assert.ok(b.equilibrium);assert.deepEqual(b.graph.vertexPositions,fixed.graph.vertexPositions.map(Math.fround));
  assert.ok(b.axialForces.every(f=>f<0));
  const incompatible={...fixed,graph:{...fixed.graph,vertexPositions:[-1,0,0,.3,.1,.6,1,0,0]}};
  ok(typedCompute(other,incompatible));b=typedResult(other);assert.equal(b.status,'reciprocation_incomplete');assert.equal(b.equilibrium,false);assert.ok(b.maxResidual>.01);
  assert.deepEqual(b.graph.vertexPositions,incompatible.graph.vertexPositions.map(Math.fround));
  for(const filename of ['force-polyhedra-chain.json','force-polyhedra-branch.json']){
    const example=JSON.parse(await readFile(new URL('../examples/'+filename,import.meta.url),'utf8'));
    ok(typedCompute(other,example));const solved=typedResult(other);assert.ok(solved.equilibrium,filename);
    for(const id of example.supports)assert.deepEqual(solved.graph.vertexPositions.slice(3*id,3*id+3),example.graph.vertexPositions.slice(3*id,3*id+3).map(Math.fround));
  }
  const foreign=call('zspace_vault_create');assert.equal(call('zspace_forcepolyhedra_result',[foreign]),0);assert.match(error(),/handle/);ok(call('zspace_toolsets_destroy',[foreign]));
  ok(compute(other,{graph:{vertexPositions:[0,0,0],edgeConnects:[]}}));assert.equal(result(other).cells[0].status,'no_forces');
}finally{ok(call('zspace_forcepolyhedra_destroy',[h]));ok(call('zspace_forcepolyhedra_destroy',[other]));}
assert.equal(call('zspace_forcepolyhedra_result',[h]),0);
console.log('Force Polyhedra WASM: EGI/AP, weighted 3DGS, cell mapping, virtual faces, force-fixed solve, typed parity, rollback and session isolation passed.');
