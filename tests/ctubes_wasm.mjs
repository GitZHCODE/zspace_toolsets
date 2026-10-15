// Run after wasm/scripts/build-wasm.bat: node tests/ctubes_wasm.mjs
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
const factory = (await import(new URL('../wasm/out/zspace_toolsets.js', import.meta.url))).default;
const module = await factory({ wasmBinary: await readFile(new URL('../wasm/out/zspace_toolsets.wasm', import.meta.url)) });
const call = (name, args = [], types = args.map(() => 'number')) => module.ccall(name, 'number', types, args);
const error = () => module.UTF8ToString(call('zspace_toolsets_last_error'));
const ok = result => assert.equal(result, 1, error());
assert.equal(call('zspace_toolsets_abi_version'), 1);
const handles = [call('zspace_ctubes_create'), call('zspace_ctubes_create')];
assert.notEqual(handles[0], handles[1]);
const input = {
  directrix: { vertexPositions: [0,0,0, 0,0,2, 1,0,4, 2,1,6], edgeConnects: [0,1,1,2,2,3] },
  generatrix: { vertexPositions: [-.25,-.25,0, .25,-.25,0, .25,.25,0, -.25,.25,0], edgeConnects: [0,1,1,2,2,3,3,0] },
  samples: 32, spline: true, mode: 0
};
for (let mode = 0; mode < 2; ++mode) {
  const h = handles[mode];
  module.FS.writeFile('/test.json', JSON.stringify({ ...input, mode }));
  ok(call('zspace_ctubes_compute', [h, '/test.json'], ['number', 'string']));
  ok(call('zspace_ctubes_select_mesh', [h, 0, 0]));
  assert.equal(call('zspace_toolsets_positions_count', [h]), 32 * 4 * 3);
  assert.equal(call('zspace_toolsets_face_count', [h]), 31 * 4);
  ok(call('zspace_ctubes_unroll', [h, .2]));
  assert.equal(call('zspace_ctubes_strip_count', [h]), 4);
  for (let i = 0; i < 4; ++i) {
    ok(call('zspace_ctubes_select_mesh', [h, 2, i]));
    const start = call('zspace_toolsets_positions', [h]) / 4;
    const count = call('zspace_toolsets_positions_count', [h]);
    const flat = module.HEAPF32.slice(start, start + count);
    assert.ok(flat.every(Number.isFinite));
    for (let j = 2; j < count; j += 3) assert.equal(flat[j], 0);
  }
}
// A failed compute must preserve the last successful session result.
module.FS.writeFile('/bad.json', JSON.stringify({ ...input, mode: 9 }));
assert.equal(call('zspace_ctubes_compute', [handles[0], '/bad.json'], ['number', 'string']), 0);
assert.match(error(), /mode/);
ok(call('zspace_ctubes_select_mesh', [handles[0], 0, 0]));
assert.equal(call('zspace_toolsets_face_count', [handles[0]]), 124);
assert.equal(call('zspace_ctubes_select_mesh', [handles[0], 2, 99]), 0);
for (const h of handles) ok(call('zspace_toolsets_destroy', [h]));
assert.equal(call('zspace_ctubes_unroll', [handles[0], .2]), 0);
assert.match(error(), /handle/);
console.log('Toolsets WASM: both C-Tube modes, unrolling, isolated sessions and error handling passed.');
// Target sections constrain ONE generic tube; apex scales are actual optimisation variables.
const fitHandle=call('zspace_ctubes_create');
const section=(z,scale)=>({vertexPositions:[-scale,-scale,z,scale,-scale,z,scale,scale,z,-scale,scale,z],edgeConnects:[0,1,1,2,2,3,3,0]});
const fitInput={directrix:{vertexPositions:[0,0,0,0,0,2,0,0,4],edgeConnects:[0,1,1,2]},generatrices:[section(4,2),section(1.3,1),section(2,1.2)],spline:false,mode:0,fitPlaneNormals:false,fitTolerance:.001};
module.FS.writeFile('/fit.json',JSON.stringify(fitInput));
ok(call('zspace_ctubes_compute',[fitHandle,'/fit.json'],['number','string']));
let report=JSON.parse(module.UTF8ToString(call('zspace_ctubes_diagnostics',[fitHandle])));
assert.equal(report.fit.active,true);assert.equal(report.fit.converged,true);
assert.equal(report.fit.solver,'L-BFGS-B');assert.equal(report.fit.gradientMethod,'Analytic chain rule');
assert.ok(Number.isFinite(report.fit.gradientNorm));assert.ok(report.fit.evaluations>0);assert.ok(report.fit.termination);
assert.deepEqual(report.fit.stations,[3,1,2]);assert.equal(report.frames.length,4);
assert.ok(report.fit.rms<report.fit.initialRms*.01);assert.ok(report.fit.apexScales.some(a=>a>1.1));
ok(call('zspace_ctubes_select_mesh',[fitHandle,0,0]));assert.equal(call('zspace_toolsets_face_count',[fitHandle]),12);
ok(call('zspace_ctubes_unroll',[fitHandle,.2]));
report=JSON.parse(module.UTF8ToString(call('zspace_ctubes_diagnostics',[fitHandle])));
assert.equal(report.stripAreas.length,4);assert.ok(report.stripAreas.every(a=>a.sourceArea>0 && a.deviationPercent<.001));
module.FS.writeFile('/fit.json',JSON.stringify({...fitInput,mode:1}));
assert.equal(call('zspace_ctubes_compute',[fitHandle,'/fit.json'],['number','string']),0);assert.match(error(),/Generic/);
ok(call('zspace_toolsets_destroy',[fitHandle]));
console.log('C-Tube target-section fitting: one tube, inserted stations, optimised apex values, frames and area metrics passed.');
const weightedHandle=call('zspace_ctubes_create');
const computeWeighted=data=>{
  module.FS.writeFile('/weighted.json',JSON.stringify(data));
  ok(call('zspace_ctubes_compute',[weightedHandle,'/weighted.json'],['number','string']));
  return JSON.parse(module.UTF8ToString(call('zspace_ctubes_diagnostics',[weightedHandle])));
};
const weightedInput={...fitInput,directrix:{vertexPositions:[0,0,0,0,0,1,0,0,2,0,0,3,0,0,4],edgeConnects:[0,1,1,2,2,3,3,4]},generatrices:[section(0,1),section(1,2),section(2,.7),section(4,2)],fitIterations:100};
const exact=computeWeighted(weightedInput).fit;
const smooth=computeWeighted({...weightedInput,apexSmoothnessWeight:10}).fit;
const roughness=a=>a.slice(1,-1).reduce((sum,v,i)=>sum+(a[i]-2*v+a[i+2])**2,0);
assert.ok(roughness(smooth.apexScales)<roughness(exact.apexScales)*.5);
assert.ok(smooth.objective<smooth.initialObjective);
const prioritised=computeWeighted({...weightedInput,apexSmoothnessWeight:10,sectionWeights:[1,1000,1,1]}).fit;
assert.ok(prioritised.sectionRms[1]<smooth.sectionRms[1]);
const anchorInput={...input,spline:false,preserveOrientation:true,generatrix:{...section(0,1),originVertex:2}};
for(const mode of [0,1]) {
  computeWeighted({...anchorInput,mode});ok(call('zspace_ctubes_select_mesh',[weightedHandle,0,0]));
  const offset=call('zspace_toolsets_positions',[weightedHandle])/4;
  assert.ok([...module.HEAPF32.slice(offset+6,offset+9)].every(v=>Math.abs(v)<1e-6));
}
const anchored=computeWeighted({...fitInput,generatrices:[{...section(0,1),originVertex:0},{...section(4,2),originVertex:0}]}).fit;
assert.ok(anchored.maxError<.001);
for(const invalid of [{targetPlaneWeight:-1},{sectionWeights:[1]},{generatrices:[{...section(0,1),originVertex:20},section(4,2)]}]) {
  module.FS.writeFile('/weighted.json',JSON.stringify({...fitInput,...invalid}));
  assert.equal(call('zspace_ctubes_compute',[weightedHandle,'/weighted.json'],['number','string']),0);
}
ok(call('zspace_toolsets_destroy',[weightedHandle]));
console.log('Weighted objectives, smoothness tradeoffs, per-target priorities and geometric anchors passed.');


// Open sections generate boundary strips without a synthetic closing edge.
const openHandle=call('zspace_ctubes_create');
const line=(z,width)=>({vertexPositions:[-width,0,z,width,0,z],edgeConnects:[0,1]});
const openInput={directrix:{vertexPositions:[0,0,0,0,0,2,0,0,4],edgeConnects:[0,1,1,2]},generatrices:[line(0,1),{vertexPositions:[2,0,4,-2,0,4],edgeConnects:[0,1]}],spline:false,fitPlaneNormals:false};
module.FS.writeFile('/open.json',JSON.stringify(openInput));ok(call('zspace_ctubes_compute',[openHandle,'/open.json'],['number','string']));
const openFit=JSON.parse(module.UTF8ToString(call('zspace_ctubes_diagnostics',[openHandle]))).fit;assert.equal(openFit.converged,true);
ok(call('zspace_ctubes_select_mesh',[openHandle,0,0]));assert.equal(call('zspace_toolsets_face_count',[openHandle]),2);
ok(call('zspace_ctubes_unroll',[openHandle,.2]));const openReport=JSON.parse(module.UTF8ToString(call('zspace_ctubes_diagnostics',[openHandle])));assert.equal(openReport.stripAreas.length,1);assert.ok(openReport.stripAreas[0].deviationPercent<.001);
ok(call('zspace_toolsets_destroy',[openHandle]));console.log('Open line targets, reversed endpoint matching, seam-free topology and unrolling passed.');

const tHandle=call('zspace_ctubes_create');
const tSection={vertexPositions:[0,0,0,-1,0,0,1,0,0,0,-1,0],edgeConnects:[0,1,0,2,0,3],originVertex:0};
for(const mode of [0,1]){
 module.FS.writeFile('/t.json',JSON.stringify({...openInput,generatrices:[tSection],mode,preserveOrientation:true}));
 ok(call('zspace_ctubes_compute',[tHandle,'/t.json'],['number','string']));
 ok(call('zspace_ctubes_select_mesh',[tHandle,0,0]));assert.equal(call('zspace_toolsets_face_count',[tHandle]),6);
 ok(call('zspace_ctubes_unroll',[tHandle,.2]));const report=JSON.parse(module.UTF8ToString(call('zspace_ctubes_diagnostics',[tHandle])));
 assert.equal(report.stripAreas.length,2);assert.ok(report.stripAreas.every(a=>a.deviationPercent<.001));
}
module.FS.writeFile('/t.json',JSON.stringify({...openInput,generatrices:[tSection,{vertexPositions:[0,-2,4,2,0,4,0,0,4,-2,0,4],edgeConnects:[2,3,1,2,0,2],originVertex:2}]}));
ok(call('zspace_ctubes_compute',[tHandle,'/t.json'],['number','string']));
assert.equal(JSON.parse(module.UTF8ToString(call('zspace_ctubes_diagnostics',[tHandle]))).fit.converged,true);
ok(call('zspace_toolsets_destroy',[tHandle]));console.log('Branched T sections, shared anchors, topology matching and continuous flange and stem unrolling passed.');

// Near-endpoint target planes formerly inserted a sliver with a collapsed stem edge.
const nearHandle=call('zspace_ctubes_create');
module.FS.writeFile('/near.json',await readFile(new URL('./fixtures/ctubes-near-endpoint.json',import.meta.url),'utf8'));
ok(call('zspace_ctubes_compute',[nearHandle,'/near.json'],['number','string']));
const nearFit=JSON.parse(module.UTF8ToString(call('zspace_ctubes_diagnostics',[nearHandle])));
assert.equal(nearFit.frames.length,65);
ok(call('zspace_ctubes_unroll',[nearHandle,.2]));
const nearAreas=JSON.parse(module.UTF8ToString(call('zspace_ctubes_diagnostics',[nearHandle]))).stripAreas;
assert.equal(nearAreas.length,2);assert.ok(nearAreas.every(a=>a.deviationPercent<.01));
ok(call('zspace_toolsets_destroy',[nearHandle]));console.log('Near-endpoint Rhino T targets with sketch smoothness weights unroll without collapsed edges.');
