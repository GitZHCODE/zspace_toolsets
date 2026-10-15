import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
const module=await(await import('../wasm/out/zspace_toolsets.js')).default({wasmBinary:await readFile(new URL('../wasm/out/zspace_toolsets.wasm',import.meta.url))});
const call=(name,args=[],types=args.map(()=> 'number'))=>module.ccall(name,'number',types,args);
const ok=value=>assert.equal(value,1,module.UTF8ToString(call('zspace_toolsets_last_error')));
const fixture=JSON.parse(await readFile(new URL('./fixtures/ctubes_branch_closed.json',import.meta.url),'utf8'));
for(let first=0;first<fixture.length;first+=3){
  const original=fixture.slice(first,first+3),h=call('zspace_ctubes_create');
  module.FS.writeFile('/central.json',JSON.stringify({branches:original.map(i=>({frames:i.frames,centralTargets:i.flangeTargets})),fixedDirectrices:true,iterations:100,outwardLength:.2,stemLength:.4,tolerance:1e-5}));
  ok(call('zspace_ctubes_branch_optimise',[h,'/central.json'],['number','string']));
  const output=JSON.parse(module.UTF8ToString(call('zspace_ctubes_branch_result',[h])));
  console.log(output.report);
  assert.ok(output.report.gradientCheckError<1e-5);
  assert.ok(output.report.converged,`Seam error ${output.report.seamError}`);
  assert.ok(output.report.seamError<=1e-5);
  const shared=new Map();
  for(const [b,r]of output.branches.entries()){
    for(let j=0;j<r.frames.length;j++){
      assert.deepEqual(r.frames[j].origin,original[b].frames[j].origin.map(Math.fround));
      const key=JSON.stringify(original[b].flangeTargets[j]),tip=r.sections.slice(12*j+6,12*j+9);
      if(shared.has(key))assert.ok(Math.hypot(...tip.map((x,k)=>x-shared.get(key)[k]))<1.1e-5);else shared.set(key,tip);
    }
    const tube=call('zspace_ctubes_create');
    module.FS.writeFile('/tube.json',JSON.stringify({...original[b],frames:r.frames,sections:r.sections,generatrix:{vertexPositions:r.sections.slice(0,12),edgeConnects:[1,0,0,2,0,3]}}));
    ok(call('zspace_ctubes_compute',[tube,'/tube.json'],['number','string']));ok(call('zspace_ctubes_unroll',[tube,.2]));assert.equal(call('zspace_ctubes_strip_count',[tube]),2);ok(call('zspace_toolsets_destroy',[tube]));
  }
  ok(call('zspace_toolsets_destroy',[h]));
  console.log('Fixed directrices: shared central seam fitted; six strips unroll.');
}
