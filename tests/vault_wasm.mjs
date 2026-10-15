import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
const factory = (await import('../wasm/out/zspace_toolsets.js')).default;
const m = await factory({ wasmBinary: await readFile(new URL('../wasm/out/zspace_toolsets.wasm', import.meta.url)) });
const call = (name,args=[],types=args.map(()=>'number')) => m.ccall(name,'number',types,args);
const error = () => m.UTF8ToString(call('zspace_toolsets_last_error'));
const ok = n => assert.equal(n,1,error());
const result = h => { const p=call('zspace_vault_result',[h]);assert.ok(p,error());return JSON.parse(m.UTF8ToString(p)); };
const fit = (h,input) => { m.FS.writeFile('/vault-test.json',JSON.stringify(input));return call('zspace_vault_fit',[h,'/vault-test.json'],['number','string']); };
const h=call('zspace_vault_create'), other=call('zspace_vault_create'), tubes=call('zspace_ctubes_create');
const input={mesh:{vertexPositions:[-1,0,0,1,0,0,0,0,1],faceCounts:[3],faceConnects:[0,1,2]},supports:[0,1],masses:[0,0,2],gravity:[0,0,-1]};
ok(fit(h,input)); let r=result(h);assert.ok(r.converged);assert.equal(r.activeEdges,2);
for(const e of r.activeEdgeIds)assert.ok(Math.abs(r.forceDensities[e]+1)<1e-9);
assert.ok(Math.abs(r.supportReactions[2]-1)<1e-6);
assert.equal(call('zspace_vault_step',[tubes,1]),0);assert.match(error(),/Vault/);
assert.equal(call('zspace_ctubes_unroll',[h,1]),0);assert.match(error(),/C-Tubes/);
ok(fit(other,{...input,settings:{qMin:.1,qMax:.2,relativeTolerance:.001}}));
const before=result(other);ok(call('zspace_vault_step',[other,10]));r=result(other);
assert.ok(r.iteration>0 && r.residualNorm<before.residualNorm);
for(const e of r.activeEdgeIds)assert.ok(r.forceDensities[e]>=-.2 && r.forceDensities[e]<=-.1);
ok(call('zspace_vault_select_mesh',[other]));const start=call('zspace_toolsets_positions',[other])/4;
assert.deepEqual(Array.from(m.HEAPF32.slice(start,start+6)),[-1,0,0,1,0,0]);
assert.equal(result(h).iteration,0); // independent sessions
ok(call('zspace_vault_previous',[other]));assert.equal(result(other).iteration,r.iteration-1);
ok(call('zspace_vault_reset',[other]));assert.deepEqual(result(other).forceDensities,before.forceDensities);
assert.equal(fit(other,{...input,masses:[-1,0,2]}),0);assert.match(error(),/mass/i);assert.equal(result(other).initialResidual,before.initialResidual);
assert.equal(fit(other,{...input,settings:{qMin:2,qMax:1}}),0);assert.match(error(),/qMin/);
assert.equal(fit(other,{...input,mesh:{...input.mesh,faceConnects:[0,1,99]}}),0);assert.match(error(),/index/);
assert.equal(fit(other,{...input,tensionEdges:[[0,99]]}),0);assert.match(error(),/edge/);
ok(fit(other,{...input,mesh:{...input.mesh,vertexPositions:[-1,0,0,1,0,0,0,0,-1]},tensionEdges:[[2,0],[1,2]]}));
r=result(other);assert.ok(r.converged);for(const e of r.activeEdgeIds)assert.ok(Math.abs(r.forceDensities[e]-1)<1e-9);
const mixed={mesh:{vertexPositions:[-1,0,-1,0,-1,1,1,0,-1,0,1,1,0,0,0],faceCounts:[3,3,3,3],faceConnects:[0,1,4,1,2,4,2,3,4,3,0,4]},
  supports:[0,1,2,3],masses:[0,0,0,0,4],gravity:[0,0,-1],tensionEdges:[[1,4],[3,4]]};
for(const geometryWeight of [.5,1]){
  const cancellation={...mixed,mesh:{...mixed.mesh,vertexPositions:[-1,0,0,0,-1,0,1,0,0,0,3,0,0,0,0]},
    masses:[0,0,0,0,1],tensionEdges:[[3,4]],settings:{qMin:1,qMax:1,geometryWeight,fixXY:true}};
  ok(fit(other,cancellation));const initial=result(other);
  ok(call('zspace_vault_step',[other,1]));const rescued=result(other);
  assert.equal(rescued.iteration,1);assert.ok(rescued.residualNorm<initial.residualNorm);
  for(const edge of rescued.activeEdgeIds)assert.equal(Math.abs(rescued.forceDensities[edge]),1);
  ok(call('zspace_vault_select_mesh',[other]));const base=call('zspace_toolsets_positions',[other])/4;
  assert.deepEqual(Array.from(m.HEAPF32.slice(base,base+14)),cancellation.mesh.vertexPositions.slice(0,14));
  ok(call('zspace_vault_previous',[other]));assert.equal(result(other).iteration,0);
  ok(call('zspace_vault_step',[other,1]));ok(call('zspace_vault_reset',[other]));
  assert.equal(result(other).residualNorm,initial.residualNorm);
}
ok(fit(other,mixed));const saved=result(other);
const dp=call('zspace_vault_diagrams',[other,1e-5]);assert.ok(dp,error());const diagrams=JSON.parse(m.UTF8ToString(dp));
assert.equal(diagrams.status,'reciprocal');assert.ok(diagrams.maxEdgeMismatch<1e-5);
assert.deepEqual(result(other).forceDensities,saved.forceDensities);
assert.equal(diagrams.force.densities.filter(q=>q>0).length,2);assert.equal(diagrams.force.densities.filter(q=>q<0).length,2);
assert.equal(diagrams.form.densities.filter(q=>q===0).length,4);
for(let i=0;i<diagrams.formToForce.length;i++){const j=diagrams.formToForce[i];if(j>=0)assert.equal(diagrams.forceToForm[j],i);}
assert.equal(call('zspace_vault_diagrams',[tubes,1e-5]),0);assert.match(error(),/Vault/);
assert.equal(call('zspace_vault_diagrams',[other,-1]),0);assert.match(error(),/tolerance/);
const ring={mesh:{vertexPositions:[-2,-2,0,2,-2,0,2,2,0,-2,2,0,-1,-1,1,1,-1,1,1,1,1,-1,1,1],faceCounts:[4,4,4,4],faceConnects:[0,1,5,4,1,2,6,5,2,3,7,6,3,0,4,7]},supports:[0,1,2,3],gravity:[0,0,-1]};
ok(fit(other,ring));const ringQ=result(other).forceDensities;
const ringPointer=call('zspace_vault_diagrams',[other,1e-5]);assert.ok(ringPointer,error());
const ringDiagrams=JSON.parse(m.UTF8ToString(ringPointer));assert.equal(ringDiagrams.exact,true);
assert.equal(ringDiagrams.force.vertexPositions.length,15);assert.equal(ringDiagrams.force.edges.length,16);
assert.deepEqual(result(other).forceDensities,ringQ);
for(const fixXY of [false,true]){
  const positions=[-1,-1,0,1,1,0,.4,.2,1];
  ok(fit(other,{...input,mesh:{...input.mesh,vertexPositions:positions},settings:{fixXY,qMin:.1,qMax:.2,relativeTolerance:1e-6}}));
  const initial=result(other);ok(call('zspace_vault_step',[other,10]));assert.ok(result(other).residualNorm<initial.residualNorm);
  ok(call('zspace_vault_select_mesh',[other]));const p=call('zspace_toolsets_positions',[other])/4;
  const xyz=Array.from(m.HEAPF32.slice(p,p+9));assert.deepEqual(xyz.slice(0,6),positions.slice(0,6));assert.notEqual(xyz[8],1);
  if(fixXY)assert.deepEqual(xyz.slice(6,8),Array.from(new Float32Array([.4,.2])));
  else {assert.notEqual(xyz[6],Math.fround(.4));assert.notEqual(xyz[7],Math.fround(.2));}
}
for(const fixXY of [false,true]){
  const reg={...input,mesh:{...input.mesh,vertexPositions:[-1,0,0,1,0,0,.6,0,1]},settings:{fixXY,minMaxEdgeRatio:.9,edgeRegularizationStrength:1}};
  ok(fit(other,reg));const initial=result(other);assert.equal(initial.converged,false);assert.ok(initial.residualNorm<1e-6);
  ok(call('zspace_vault_step',[other,20]));const after=result(other);
  assert.ok(after.minMaxEdgeRatio>initial.minMaxEdgeRatio);assert.ok(after.edgeRegularizationError<initial.edgeRegularizationError);
  ok(call('zspace_vault_select_mesh',[other]));const p=call('zspace_toolsets_positions',[other])/4;
  assert.deepEqual(Array.from(m.HEAPF32.slice(p,p+6)),[-1,0,0,1,0,0]);
  if(fixXY)assert.equal(m.HEAPF32[p+6],Math.fround(.6));
  assert.equal(fit(other,{...reg,settings:{minMaxEdgeRatio:1.1}}),0);assert.match(error(),/ratio/);
}
let freeDrift=0;
for(const planConstraint of [0,1,2]){
  const plan={...input,mesh:{...input.mesh,vertexPositions:[-1,0,0,1,0,0,.6,0,1]},settings:{planConstraint,xyPenaltyWeight:20,minMaxEdgeRatio:.9,edgeRegularizationStrength:1}};
  ok(fit(other,plan));ok(call('zspace_vault_step',[other,20]));const r=result(other);
  if(planConstraint===0){freeDrift=r.rmsXYDisplacement;assert.ok(freeDrift>0);}
  if(planConstraint===1)assert.ok(r.rmsXYDisplacement>0&&r.rmsXYDisplacement<freeDrift);
  if(planConstraint===2)assert.equal(r.rmsXYDisplacement,0);
}
const refInput={...input,mesh:{...input.mesh,vertexPositions:[-1,0,0,1,0,0,.6,0,1]},referencePositions:[-1,0,0,1,0,0,.2,0,1],settings:{planConstraint:1,xyPenaltyWeight:2}};
ok(fit(other,refInput));const refStart=result(other);assert.ok(refStart.rmsXYDisplacement>.39);assert.equal(refStart.converged,false);
ok(call('zspace_vault_step',[other,10]));assert.ok(result(other).rmsXYDisplacement<refStart.rmsXYDisplacement);
ok(call('zspace_vault_reset',[other]));assert.equal(result(other).rmsXYDisplacement,refStart.rmsXYDisplacement);
assert.equal(fit(other,{...refInput,referencePositions:[1,2,3]}),0);assert.match(error(),/Reference/);
for(const handle of [h,other,tubes])ok(call('zspace_toolsets_destroy',[handle]));
assert.equal(call('zspace_vault_step',[h,1]),0);assert.match(error(),/handle/);
console.log('Vault WASM: analytic native parity, masses, signs, bounds, history, session isolation and atomic failed fit passed.');
