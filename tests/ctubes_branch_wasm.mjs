import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
const module=await (await import('../wasm/out/zspace_toolsets.js')).default({wasmBinary:await readFile(new URL('../wasm/out/zspace_toolsets.wasm',import.meta.url))});
const call=(name,args=[],types=args.map(()=> 'number'))=>module.ccall(name,'number',types,args);
const error=()=>module.UTF8ToString(call('zspace_toolsets_last_error'));
const ok=value=>assert.equal(value,1,error());
const inputs=[...JSON.parse(await readFile(new URL('./fixtures/ctubes_branch_frames.json',import.meta.url),'utf8')),...JSON.parse(await readFile(new URL('./fixtures/ctubes_branch_closed.json',import.meta.url),'utf8'))];
for(const [i,input] of inputs.entries()){
  const h=call('zspace_ctubes_create');module.FS.writeFile('/branch.json',JSON.stringify(input));
  ok(call('zspace_ctubes_compute',[h,'/branch.json'],['number','string']));
  const report=JSON.parse(module.UTF8ToString(call('zspace_ctubes_diagnostics',[h])));
  assert.equal(report.frames.length,input.frames.length);assert.equal(report.fit.active,false);
  ok(call('zspace_ctubes_select_mesh',[h,0,0]));const offset=call('zspace_toolsets_positions',[h])/4,count=call('zspace_toolsets_positions_count',[h]);const p=module.HEAPF32.slice(offset,offset+count);
  assert.ok(p.every(Number.isFinite));
  // Branched display mesh duplicates the endpoints of each seed edge per station.
  for(let j=0;j<input.frames.length;j++)for(let k=0;k<6;k++){
    const f=input.frames[j],index=2*input.frames.length*Math.floor(k/2)+2*j+k%2,v=Array.from(p.slice(3*index,3*index+3));
    const distance=v.reduce((s,x,a)=>s+(x-f.origin[a])*f.normal[a],0);assert.ok(Math.abs(distance)<1e-4,`branch ${i} station ${j}: ${distance}`);
  }
  if(input.flangeTargets)for(let j=0;j<input.frames.length;j++){const index=2*input.frames.length+2*j+1,tip=Array.from(p.slice(3*index,3*index+3));assert.ok(Math.hypot(...tip.map((v,k)=>v-input.flangeTargets[j][k]))<1e-6,'Inward flange must reach its central target.');}
  const unfolded=call('zspace_ctubes_unroll',[h,.2]);
  if(!unfolded&&input.flangeTargets){assert.match(error(),/without a seam or metric distortion/);console.log(`Closed boundary ${i+1}: exact targets satisfied; strict flattening needs a seam or distortion.`);}
  else ok(unfolded);
  assert.equal(call('zspace_ctubes_strip_count',[h]),2);
  const areas=JSON.parse(module.UTF8ToString(call('zspace_ctubes_diagnostics',[h]))).stripAreas;
  if(unfolded){assert.equal(areas.length,2);assert.ok(areas.every(a=>a.sourceArea>0&&Number.isFinite(a.deviationPercent)));}else assert.equal(areas.length,0);
  console.log(`Boundary ${i+1}: ${input.frames.length} fixed planes, 2 strips, area deviations ${areas.map(a=>a.deviationPercent.toFixed(5)).join(', ')}%.`);
  module.FS.writeFile('/bad-frame.json',JSON.stringify({...input,frames:input.frames.slice(1)}));assert.equal(call('zspace_ctubes_compute',[h,'/bad-frame.json'],['number','string']),0);assert.match(error(),/one frame/);
  ok(call('zspace_ctubes_select_mesh',[h,0,0]));ok(call('zspace_toolsets_destroy',[h]));
}
// Fit boundary directrices with paper-style conical projection while retaining
// the identical shared central seams. Test the actual saved Y geometry.
for(let first=3;first<inputs.length;first+=3){
  const originals=inputs.slice(first,first+3),h=call('zspace_ctubes_create');
  module.FS.writeFile('/optimise.json',JSON.stringify({branches:originals.map(i=>({frames:i.frames,centralTargets:i.flangeTargets})),iterations:100,outwardLength:.2,stemLength:.4,tolerance:.001}));
  ok(call('zspace_ctubes_branch_optimise',[h,'/optimise.json'],['number','string']));
  const result=JSON.parse(module.UTF8ToString(call('zspace_ctubes_branch_result',[h])));
  assert.ok(result.report.gradientCheckError<1e-5);
  assert.ok(result.report.objective<result.report.initialObjective);
  for(const [b,r] of result.branches.entries()){
    const original=originals[b],tube=call('zspace_ctubes_create');
    for(let j=0;j<r.frames.length;j++)assert.ok(Math.hypot(...r.sections.slice(12*j+6,12*j+9).map((x,k)=>x-original.flangeTargets[j][k]))<1e-6);
    module.FS.writeFile('/optimised.json',JSON.stringify({...original,directrix:{...original.directrix,vertexPositions:r.frames.flatMap(f=>f.origin)},generatrix:{vertexPositions:r.sections.slice(0,12),edgeConnects:[1,0,0,2,0,3],originVertex:0},frames:r.frames,sections:r.sections}));
    ok(call('zspace_ctubes_compute',[tube,'/optimised.json'],['number','string']));
    ok(call('zspace_ctubes_unroll',[tube,.2]));
    assert.equal(call('zspace_ctubes_strip_count',[tube]),2);
    const areas=JSON.parse(module.UTF8ToString(call('zspace_ctubes_diagnostics',[tube]))).stripAreas;
    assert.equal(areas.length,2);assert.ok(areas.every(a=>a.deviationPercent<.001));
    ok(call('zspace_toolsets_destroy',[tube]));
  }
  console.log(`Optimised Y: RMS movement ${result.report.rmsDisplacement}; gradient error ${result.report.gradientCheckError}; all six strips unroll.`);
  ok(call('zspace_toolsets_destroy',[h]));
}
