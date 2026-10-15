import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
const m=await(await import('../wasm/out/zspace_toolsets.js')).default({wasmBinary:await readFile(new URL('../wasm/out/zspace_toolsets.wasm',import.meta.url))});
const call=(name,args=[],types=args.map(()=> 'number'))=>m.ccall(name,'number',types,args);
const ok=value=>assert.equal(value,1,m.UTF8ToString(call('zspace_toolsets_last_error')));
for(const fixedDirectrices of [true,false]){
  // One lifted station creates only two nonplanar panels on each ribbon.
  const branches=[-1,1,2].map(x=>({frames:Array.from({length:7},(_,j)=>({origin:[x,j,!fixedDirectrices&&j===3?.04:0],x:[-Math.sign(x),0,0],y:[0,0,1],normal:[0,Math.sign(x),0]})),centralTargets:Array.from({length:7},(_,j)=>[0,j,fixedDirectrices&&j===3?.04:0])}));
  const h=call('zspace_ctubes_create');
  m.FS.writeFile('/fit.json',JSON.stringify({branches,fixedDirectrices,fixXY:true,iterations:100,outwardLength:.2,stemLength:.4,tolerance:1e-5}));
  ok(call('zspace_ctubes_branch_optimise',[h,'/fit.json'],['number','string']));
  const output=JSON.parse(m.UTF8ToString(call('zspace_ctubes_branch_result',[h])));
  assert.ok(output.report.converged);assert.ok(output.report.gradientCheckError<1e-5);
  assert.ok(output.report.maxDisplacement<.041,'A small local defect should not substantially reshape the graph');
  for(const [b,r] of output.branches.entries()){
    for(let j=0;j<r.frames.length;j++){
      const before=branches[b],tip=r.sections.slice(12*j+6,12*j+9);
      assert.deepEqual(r.frames[j].origin.slice(0,2),before.frames[j].origin.slice(0,2));
      assert.deepEqual(tip.slice(0,2),before.centralTargets[j].slice(0,2));
      if(fixedDirectrices)assert.deepEqual(r.frames[j].origin,before.frames[j].origin);
      else assert.deepEqual(tip,before.centralTargets[j]);
    }
    m.FS.writeFile('/tube.json',JSON.stringify({directrix:{vertexPositions:r.frames.flatMap(f=>f.origin),edgeConnects:Array.from({length:6},(_,j)=>[j,j+1]).flat()},generatrix:{vertexPositions:r.sections.slice(0,12),edgeConnects:[1,0,0,2,0,3]},frames:r.frames,sections:r.sections,spline:false,mode:0}));
    ok(call('zspace_ctubes_compute',[h,'/tube.json'],['number','string']));ok(call('zspace_ctubes_unroll',[h,.2]));assert.equal(call('zspace_ctubes_strip_count',[h]),2);
  }
  console.log(fixedDirectrices?'Central graph':'Boundary directrices',output.report);
  ok(call('zspace_toolsets_destroy',[h]));
}
