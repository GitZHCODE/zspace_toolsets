import {readFile,writeFile} from 'node:fs/promises';
import {parseOBJ,extractBranchGuides,spacedBranchGuides,tProfile} from '../../../zspace_alice_webviewer/sketches/ctubeBranchInput.js';
const base=extractBranchGuides(parseOBJ(await readFile(new URL('../../../zspace_alice_webviewer/tests/fixtures/CTube_01_Y.obj',import.meta.url),'utf8')));
const fixture=[.4,.65].flatMap(spacing=>spacedBranchGuides(base,spacing).directrices.map(d=>({directrix:{vertexPositions:d.vertexPositions,edgeConnects:d.edgeConnects},generatrix:tProfile(d.frames[0],.4,.4),frames:d.frames.map(({origin,x,y,normal})=>({origin,x,y,normal})),flangeTargets:d.frames.map(f=>f.target),flangeVertex:2,spline:false,mode:0,apexScale:1})));
await writeFile(new URL('./ctubes_branch_closed.json',import.meta.url),JSON.stringify(fixture));
console.log('Updated station counts',fixture.map(f=>f.frames.length));
