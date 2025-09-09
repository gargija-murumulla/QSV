// ---------- Utilities ----------
const ZERO_C = math.complex(0,0);
function c(re, im=0) { return math.complex(re, im); }
const cre = math.re, cim = math.im;

// ---------- Gate matrices (2x2) ----------
const SQRT1_2 = 1/Math.sqrt(2);
const GATES = {
  X: [[c(0,0), c(1,0)], [c(1,0), c(0,0)]],
  Y: [[c(0,0), c(0,-1)], [c(0,1), c(0,0)]],
  Z: [[c(1,0), c(0,0)], [c(0,0), c(-1,0)]],
  H: [[c(SQRT1_2,0), c(SQRT1_2,0)], [c(SQRT1_2,0), c(-SQRT1_2,0)]],
  S: [[c(1,0), c(0,0)], [c(0,0), c(0,1)]],                      // diag(1, i)
  Sdg: [[c(1,0), c(0,0)], [c(0,0), c(0,-1)]],                   // diag(1, -i)
  T: [[c(1,0), c(0,0)], [c(0,0), math.exp(math.multiply(c(0,1), Math.PI/4))]], // diag(1, e^{iπ/4})
  Tdg: [[c(1,0), c(0,0)], [c(0,0), math.exp(math.multiply(c(0,-1), Math.PI/4))]],
};

// Parameterized single-qubit gates
function Rx(theta){
  const th = theta/2;
  return [
    [c(Math.cos(th),0), c(0,-Math.sin(th))],
    [c(0,-Math.sin(th)), c(Math.cos(th),0)]
  ];
}
function Ry(theta){
  const th = theta/2;
  return [
    [c(Math.cos(th),0), c(-Math.sin(th),0)],
    [c(Math.sin(th),0), c(Math.cos(th),0)]
  ];
}
function Rz(theta){
  const th = theta/2;
  return [
    [math.exp(math.multiply(c(0,-1), th)), c(0,0)],
    [c(0,0), math.exp(math.multiply(c(0,1), th))]
  ];
}
function Phase(phi){
  return [[c(1,0), c(0,0)], [c(0,0), math.exp(math.multiply(c(0,1), phi))]];
}

// ---------- DOM elements ----------
const btnSet = document.getElementById('btnSet');
const afterSet = document.getElementById('afterSet');
const numQInput = document.getElementById('numQ');
const basisSelect = document.getElementById('basisState');

const gateType = document.getElementById('gateType');
const singleTargetDiv = document.getElementById('singleTargetDiv');
const targetQ = document.getElementById('targetQ');

const angleDiv = document.getElementById('angleDiv');
const angleDeg = document.getElementById('angleDeg');
const angleLabel = document.getElementById('angleLabel');

const cnotDiv = document.getElementById('cnotDiv');
const controlQ = document.getElementById('controlQ');
const targetQ2 = document.getElementById('targetQ2');

const swapDiv = document.getElementById('swapDiv');
const swapA = document.getElementById('swapA');
const swapB = document.getElementById('swapB');

const ccnotDiv = document.getElementById('ccnotDiv');
const cc_c1 = document.getElementById('cc_c1');
const cc_c2 = document.getElementById('cc_c2');
const cc_t = document.getElementById('cc_t');

const btnAddGate = document.getElementById('btnAddGate');
const btnUndo = document.getElementById('btnUndo');
const btnClearGates = document.getElementById('btnClearGates');
const gatesListDiv = document.getElementById('gatesList');
const btnRun = document.getElementById('btnRun');

const resultsDiv = document.getElementById('results');
const blochSpheresDiv = document.getElementById('blochSpheres');
const container = document.getElementById("quantumCircuit");
  
// ---------- App state ----------
let nQ = 2;
let stateVec = []; // array of math.complex
let gateSequence = []; // each gate object: {type, params, angle?}

// ---------- Setup handlers ----------
btnSet.addEventListener('click', onSet);
gateType.addEventListener('change', onGateTypeChange);
btnAddGate.addEventListener('click', onAddGate);
btnUndo.addEventListener('click', onUndo);
btnClearGates.addEventListener('click', onClearGates);
btnRun.addEventListener('click', onRun);

// initialize UI
onGateTypeChange();

// ---------- Functions ----------
blochSpheresDiv.innerHTML = "<div class = grid > <b><h2 style= text-align:'centre'>Qubit</h2></b> <p>The fundemental unit of quantum information, serving as the quantum equvivalent of a classical computer's bit. A qubit can have states 0, 1, 0/1(superposition). </p></div>"
    
let firstQubit = false;
function onSet(){
  nQ = parseInt(numQInput.value);
  if (!(nQ >=1 && nQ <=5)) { alert("Choose n between 1 and 5"); return; }
  populateBasis(nQ);
  populateQubitSelectors(nQ);
  initState(nQ);
  afterSet.classList.remove('hidden');
  gateSequence = [];
  renderGateList();
  resultsDiv.innerHTML = "<div class='small'>Initial state set. Add gates and click Run.</div>";
  if(!firstQubit){
    blochSpheresDiv.innerHTML = "<div class = grid ><p>Tensor  products (&#8855;) are essential for describing subsystems composed of multiple quantum subsystems, where the state of the total system is given by the tensor product of the states of the individual subsystems </p></div>";
    firstQubit = true;
  }
  blochSpheresDiv.innerHTML = "";
  resultsDiv.innerHTLML = "";
  container.innerHTML = "<h2>Circuit Diagram </h2>";
}
function populateBasis(n){
  basisSelect.innerHTML = "";
  for (let i = 0; i < (1 << n); i++){
    const opt = document.createElement('option');
    opt.value = i.toString(2).padStart(n, '0');
    // separate each qubit with | >
    opt.innerHTML = opt.value.split('').map(bit => `|${bit}⟩ `).join(' &#8855; ');
    basisSelect.appendChild(opt);
  }
  // default to |0⟩|0⟩...|0⟩
  basisSelect.value = '0'.repeat(n);
}


function populateQubitSelectors(n){
  const sels = [targetQ, controlQ, targetQ2, swapA, swapB, cc_c1, cc_c2, cc_t];
  sels.forEach(s => s.innerHTML = '');
  for (let i=0;i<n;i++){
    const opt = (id)=>{ const o=document.createElement('option'); o.value=i; o.text='q'+i; return o; };
    sels.forEach(s => s.appendChild(opt()));
  }
}

function initState(n){
  const dim = 1<<n;
  stateVec = Array(dim).fill(0).map(()=>c(0,0));
  const initIndex = parseInt(basisSelect.value || "0", 2);
  stateVec[initIndex] = c(1,0);
}

function onGateTypeChange(){
  const type = gateType.value;
  // hide all
  singleTargetDiv.classList.add('hidden');
  cnotDiv.classList.add('hidden');
  swapDiv.classList.add('hidden');
  ccnotDiv.classList.add('hidden');
  angleDiv.classList.add('hidden');

  // show relevant
  if (['X','Y','Z','H','S','Sdg','T','Tdg','Rx','Ry','Rz','Phase','MEASURE'].includes(type)){
    singleTargetDiv.classList.remove('hidden');
  }
  if (['Rx','Ry','Rz','Phase'].includes(type)){
    angleDiv.classList.remove('hidden');
    angleLabel.textContent = (type==='Phase') ? 'φ (degrees):' : 'θ (degrees):';
  }
  if (type === 'CNOT' || type === 'CZ'){
    cnotDiv.classList.remove('hidden');
  }
  if (type === 'SWAP'){
    swapDiv.classList.remove('hidden');
  }
  if (type === 'CCNOT'){
    ccnotDiv.classList.remove('hidden');
  }
}

function onAddGate(){
  const type = gateType.value;
  let gate = { type, params: [] };

  if (['X','Y','Z','H','S','Sdg','T','Tdg','Rx','Ry','Rz','Phase','MEASURE'].includes(type)){
    const t = parseInt(targetQ.value);
    gate.params = [t];
    if (['Rx','Ry','Rz','Phase'].includes(type)){
      gate.angle = (parseFloat(angleDeg.value) || 0) * Math.PI/180; // store radians
    }
  } else if (type === 'CNOT' || type === 'CZ'){
    const c = parseInt(controlQ.value), t = parseInt(targetQ2.value);
    if (c === t) { alert("Control and target must be different"); return; }
    gate.params = [c, t];
  } else if (type === 'SWAP'){
    const a = parseInt(swapA.value), b = parseInt(swapB.value);
    if (a === b) { alert("Choose two different qubits"); return; }
    gate.params = [a, b];
  } else if (type === 'CCNOT'){
    const c1 = parseInt(cc_c1.value), c2 = parseInt(cc_c2.value), t = parseInt(cc_t.value);
    const set = new Set([c1,c2,t]);
    if (set.size < 3) { alert("Controls and target must be all different"); return; }
    if (nQ < 3) { alert("CCNOT needs at least 3 qubits"); return; }
    gate.params = [c1, c2, t];
  }

  gateSequence.push(gate);
  renderGateList();
  addGate(gate);
}
function addGate(gate) {
  renderCircuit(nQ, gateSequence);
}

//circuit printing
function renderCircuit(numQubits, gates) {
  
  const numClassical = numQubits;
  container.innerHTML = "";
  container.innerHTML += "<h2>circuit Diagram </h2>";

  const width = 120 * (gates.length + 1);
  const qheight = 60 ;
  const cHeight = 40;
  const height = numQubits * qheight + numClassical * cHeight + 60;

  const svgNS = "http://www.w3.org/2000/svg";
  const svg = document.createElementNS(svgNS, "svg");
  svg.setAttribute("width", width);
  svg.setAttribute("height", height);

  // --- Draw wires ---
  for (let q = 0; q < numQubits; q++) {
    const line = document.createElementNS(svgNS, "line");
    line.setAttribute("x1", 20);
    line.setAttribute("y1", 30 + q * qheight);
    line.setAttribute("x2", width - 20);
    line.setAttribute("y2", 30 + q * qheight);
    line.setAttribute("stroke", "black");
    line.setAttribute("stroke-width", "2");
    svg.appendChild(line);

    // Label
    const text = document.createElementNS(svgNS, "text");
    text.setAttribute("x", 0);
    text.setAttribute("y", 35 + q * qheight);
    text.textContent = `q${q}`;
    svg.appendChild(text);
  }
  for (let q = 0; q < numQubits; q++) {
    const y = 30 + q * qheight;
    const line = document.createElementNS(svgNS, "line");
    line.setAttribute("x1", 20);
    line.setAttribute("y1", y);
    line.setAttribute("x2", width - 20);
    line.setAttribute("y2", y);
    line.setAttribute("stroke", "black");
    line.setAttribute("stroke-width", "2");
    svg.appendChild(line);

    // Quantum labels
    const text = document.createElementNS(svgNS, "text");
    text.setAttribute("x", 0);
    text.setAttribute("y", y + 5);
    text.textContent = `q${q}`;
    svg.appendChild(text);
  }
  const startY = numQubits * qheight + 50;
  // --- Draw classical registers ---
  for (let c = 0; c < numClassical; c++) {
    const y = startY + c * cHeight;
    const line = document.createElementNS(svgNS, "line");
    line.setAttribute("x1", 20);
    line.setAttribute("y1", y);
    line.setAttribute("x2", width - 20);
    line.setAttribute("y2", y);
    line.setAttribute("stroke", "blue");
    line.setAttribute("stroke-width", "2");
    svg.appendChild(line);

    // Classical labels
    const text = document.createElementNS(svgNS, "text");
    text.setAttribute("x", 0);
    text.setAttribute("y", y + 5);
    text.textContent = `cr[${c}]`;
    svg.appendChild(text);
  }


  // --- Draw gates ---
  gates.forEach((g, i) => {
    const x = 100 + i * 120;

    // 🎯 Handle single-qubit standard + rotation gates
    if (["X", "Y", "Z", "H", "S", "T", "Sdg", "Tdg", "Rx", "Ry", "Rz", "Phase"].includes(g.type)) {
      const y = 30 + g.params[0] * qheight;

      const rect = document.createElementNS(svgNS, "rect");
      rect.setAttribute("x", x - 25);
      rect.setAttribute("y", y - 25);
      rect.setAttribute("width", 50);
      rect.setAttribute("height", 50);
      rect.setAttribute("fill", "#d1e7dd");
      rect.setAttribute("stroke", "black");
      svg.appendChild(rect);

      const label = document.createElementNS(svgNS, "text");
      label.setAttribute("x", x);
      label.setAttribute("y", y);
      label.setAttribute("text-anchor", "middle");
      label.setAttribute("font-size", "14");
      label.setAttribute("dominant-baseline", "middle");

      // Show gate + angle if rotation
      if (["Rx", "Ry", "Rz", "Phase"].includes(g.type)) {
        const angleDeg = g.angle ? (g.angle * 180 / Math.PI).toFixed(1) : "";
        label.textContent = `${g.type}${angleDeg ? `(${angleDeg}°)` : ""}`;
        
      } else {
        label.textContent = g.type;
      }

      svg.appendChild(label);
    }

    // CNOT
    if (g.type === "CNOT") {
      const c = g.params[0];
      const t = g.params[1];
      const yc = 30 + c * qheight;
      const yt = 30 + t * qheight;

      const dot = document.createElementNS(svgNS, "circle");
      dot.setAttribute("cx", x);
      dot.setAttribute("cy", yc);
      dot.setAttribute("r", 6);
      dot.setAttribute("fill", "black");
      svg.appendChild(dot);

      const circle = document.createElementNS(svgNS, "circle");
      circle.setAttribute("cx", x);
      circle.setAttribute("cy", yt);
      circle.setAttribute("r", 12);
      circle.setAttribute("stroke", "black");
      circle.setAttribute("fill", "white");
      svg.appendChild(circle);

      const lineV = document.createElementNS(svgNS, "line");
      lineV.setAttribute("x1", x);
      lineV.setAttribute("y1", yc);
      lineV.setAttribute("x2", x);
      lineV.setAttribute("y2", yt);
      lineV.setAttribute("stroke", "black");
      lineV.setAttribute("stroke-width", "2");
      svg.appendChild(lineV);

      const lineH = document.createElementNS(svgNS, "line");
      lineH.setAttribute("x1", x - 10);
      lineH.setAttribute("y1", yt);
      lineH.setAttribute("x2", x + 10);
      lineH.setAttribute("y2", yt);
      lineH.setAttribute("stroke", "black");
      lineH.setAttribute("stroke-width", "2");
      svg.appendChild(lineH);

      const lineV2 = document.createElementNS(svgNS, "line");
      lineV2.setAttribute("x1", x);
      lineV2.setAttribute("y1", yt - 10);
      lineV2.setAttribute("x2", x);
      lineV2.setAttribute("y2", yt + 10);
      lineV2.setAttribute("stroke", "black");
      lineV2.setAttribute("stroke-width", "2");
      svg.appendChild(lineV2);
    }
    //cz
    // CZ
if (g.type === "CZ") {
  const c = g.params[0];
  const t = g.params[1];
  const yc = 30 + c * qheight;
  const yt = 30 + t * qheight;

  // Control dot
  const dotC = document.createElementNS(svgNS, "circle");
  dotC.setAttribute("cx", x);
  dotC.setAttribute("cy", yc);
  dotC.setAttribute("r", 6);
  dotC.setAttribute("fill", "black");
  svg.appendChild(dotC);

  // Target dot
  const dotT = document.createElementNS(svgNS, "circle");
  dotT.setAttribute("cx", x);
  dotT.setAttribute("cy", yt);
  dotT.setAttribute("r", 6);
  dotT.setAttribute("fill", "black");
  svg.appendChild(dotT);

  // Vertical line connecting them
  const lineV = document.createElementNS(svgNS, "line");
  lineV.setAttribute("x1", x);
  lineV.setAttribute("y1", yc);
  lineV.setAttribute("x2", x);
  lineV.setAttribute("y2", yt);
  lineV.setAttribute("stroke", "black");
  lineV.setAttribute("stroke-width", "2");
  svg.appendChild(lineV);
}

    // SWAP
    if (g.type === "SWAP") {
      const a = g.params[0];
      const b = g.params[1];
      const ya = 30 + a * qheight;
      const yb = 30 + b * qheight;

      const line1 = document.createElementNS(svgNS, "line");
      line1.setAttribute("x1", x - 10);
      line1.setAttribute("y1", ya - 10);
      line1.setAttribute("x2", x + 10);
      line1.setAttribute("y2", ya + 10);
      line1.setAttribute("stroke", "black");
      line1.setAttribute("stroke-width", "2");
      svg.appendChild(line1);

      const line2 = document.createElementNS(svgNS, "line");
      line2.setAttribute("x1", x - 10);
      line2.setAttribute("y1", ya + 10);
      line2.setAttribute("x2", x + 10);
      line2.setAttribute("y2", ya - 10);
      line2.setAttribute("stroke", "black");
      line2.setAttribute("stroke-width", "2");
      svg.appendChild(line2);

      const line3 = document.createElementNS(svgNS, "line");
      line3.setAttribute("x1", x - 10);
      line3.setAttribute("y1", yb - 10);
      line3.setAttribute("x2", x + 10);
      line3.setAttribute("y2", yb + 10);
      line3.setAttribute("stroke", "black");
      line3.setAttribute("stroke-width", "2");
      svg.appendChild(line3);

      const line4 = document.createElementNS(svgNS, "line");
      line4.setAttribute("x1", x - 10);
      line4.setAttribute("y1", yb + 10);
      line4.setAttribute("x2", x + 10);
      line4.setAttribute("y2", yb - 10);
      line4.setAttribute("stroke", "black");
      line4.setAttribute("stroke-width", "2");
      svg.appendChild(line4);

      const lineV = document.createElementNS(svgNS, "line");
      lineV.setAttribute("x1", x);
      lineV.setAttribute("y1", ya);
      lineV.setAttribute("x2", x);
      lineV.setAttribute("y2", yb);
      lineV.setAttribute("stroke", "black");
      lineV.setAttribute("stroke-width", "2");
      svg.appendChild(lineV);
    }

    // Toffoli (CCNOT)
    if (g.type === "CCNOT") {
      const c1 = g.params[0];
      const c2 = g.params[1];
      const t = g.params[2];
      const y1 = 30 + c1 * qheight;
      const y2 = 30 + c2 * qheight;
      const yt = 30 + t * qheight;

      [y1, y2].forEach(yc => {
        const dot = document.createElementNS(svgNS, "circle");
        dot.setAttribute("cx", x);
        dot.setAttribute("cy", yc);
        dot.setAttribute("r", 6);
        dot.setAttribute("fill", "black");
        svg.appendChild(dot);
      });

      const lineV = document.createElementNS(svgNS, "line");
      lineV.setAttribute("x1", x);
      lineV.setAttribute("y1", Math.min(y1, y2));
      lineV.setAttribute("x2", x);
      lineV.setAttribute("y2", yt);
      lineV.setAttribute("stroke", "black");
      lineV.setAttribute("stroke-width", "2");
      svg.appendChild(lineV);

      const circle = document.createElementNS(svgNS, "circle");
      circle.setAttribute("cx", x);
      circle.setAttribute("cy", yt);
      circle.setAttribute("r", 12);
      circle.setAttribute("stroke", "black");
      circle.setAttribute("fill", "white");
      svg.appendChild(circle);

      const lineH = document.createElementNS(svgNS, "line");
      lineH.setAttribute("x1", x - 10);
      lineH.setAttribute("y1", yt);
      lineH.setAttribute("x2", x + 10);
      lineH.setAttribute("y2", yt);
      lineH.setAttribute("stroke", "black");
      lineH.setAttribute("stroke-width", "2");
      svg.appendChild(lineH);

      const lineV2 = document.createElementNS(svgNS, "line");
      lineV2.setAttribute("x1", x);
      lineV2.setAttribute("y1", yt - 10);
      lineV2.setAttribute("x2", x);
      lineV2.setAttribute("y2", yt + 10);
      lineV2.setAttribute("stroke", "black");
      lineV2.setAttribute("stroke-width", "2");
      svg.appendChild(lineV2);
    }
    if (g.type === "MEASURE") {
      const q = g.params[0];
      const c = g.params[0];
      const yq = 30 + q * qheight;
      const yc = startY + c * cHeight;

      const rect = document.createElementNS(svgNS, "rect");
      rect.setAttribute("x", x - 20);
      rect.setAttribute("y", yq - 20);
      rect.setAttribute("width", 40);
      rect.setAttribute("height", 40);
      rect.setAttribute("fill", "#fff3cd");
      rect.setAttribute("stroke", "black");
      svg.appendChild(rect);

      const label = document.createElementNS(svgNS, "text");
      label.setAttribute("x", x);
      label.setAttribute("y", yq);
      label.setAttribute("text-anchor", "middle");
      label.setAttribute("dominant-baseline", "middle");
      label.setAttribute("font-size", "12");
      label.textContent = "M";
      svg.appendChild(label);

      const line = document.createElementNS(svgNS, "line");
      line.setAttribute("x1", x);
      line.setAttribute("y1", yq +20);
      line.setAttribute("x2", x);
      line.setAttribute("y2", yc);
      line.setAttribute("stroke", "blue");
      line.setAttribute("stroke-dasharray", "4");
      svg.appendChild(line);
    }
          // --- Draw identity gates for qubits not affected by this gate ---
    for (let q = 0; q < numQubits; q++) {
      let isTarget = false;

      if (["X","Y","Z","H","S","T","Sdg","Tdg","Rx","Ry","Rz","Phase"].includes(g.type)) {
        isTarget = (q === g.params[0]);
      } else if (["CNOT", "CZ"].includes(g.type)) {
        isTarget = (q === g.params[0] || q === g.params[1]);
      } else if (g.type === "CCNOT") {
        isTarget = (q === g.params[0] || q === g.params[1] || q === g.params[2]);
      } else if (g.type === "MEASURE") {
        isTarget = (q === g.params[0]);
      }

      if (!isTarget) {
        const y = 30 + q * qheight;
        const rect = document.createElementNS(svgNS, "rect");
        rect.setAttribute("x", x - 15);
        rect.setAttribute("y", y - 15);
        rect.setAttribute("width", 30);
        rect.setAttribute("height", 30);
        rect.setAttribute("fill", "#f0f0f0");  // gray/light color for identity
        rect.setAttribute("stroke", "black");
        svg.appendChild(rect);

        const label = document.createElementNS(svgNS, "text");
        label.setAttribute("x", x);
        label.setAttribute("y", y);
        label.setAttribute("text-anchor", "middle");
        label.setAttribute("dominant-baseline", "middle");
        label.setAttribute("font-size", "12");
        label.textContent = "I";
        svg.appendChild(label);
      }
    }
  });



  container.appendChild(svg);
}


function onUndo(){
  gateSequence.pop();
  renderGateList();
  container.innerHTML = "<h2> Circuit Diagram </h2>";
  blochSpheresDiv.innerHTML = "";
  resultsDiv.innerHTML = "";
  renderCircuit(nQ,gateSequence);
}

function onClearGates(){
  gateSequence = [];
  renderGateList();
  
  container.innerHTML = "<h2> Circuit Diagram </h2>";
  blochSpheresDiv.innerHTML = "";
  resultsDiv.innerHTML = "";
  
}

function renderGateList(){
  gatesListDiv.innerHTML = "";
  if (gateSequence.length === 0){
    gatesListDiv.innerHTML = "<div class='small'>No gates added yet.</div>";
    return;
  }
  gateSequence.forEach((g,i)=>{
    const d = document.createElement('div');
    d.className = "gate-item";

    const left = document.createElement('div');
    left.className = 'gate-left';
    let desc = `${i+1}. ${g.type}`;
    if (g.params?.length){
      desc += ` (${g.params.join(',')})`;
    }
    if (g.angle !== undefined){
      const deg = (g.angle*180/Math.PI).toFixed(2);
      desc += `, ${deg}°`;
    }
    left.textContent = desc;

    const right = document.createElement('div');
    right.className = 'gate-right';

    const up = document.createElement('button');
    up.textContent = '↑';
    up.onclick = ()=>{ if(i>0){ [gateSequence[i-1],gateSequence[i]]=[gateSequence[i],gateSequence[i-1]]; 
      renderGateList();
      container.innerHTML = "<h2>Circuit Diagram</h2>";
      renderCircuit(nQ,gateSequence);
      blochSpheresDiv.innerHTML = "";
      resultsDiv.innerHTML = "";
      } };

    const down = document.createElement('button');
    down.textContent = '↓';
    down.onclick = ()=>{ if(i<gateSequence.length-1){ [gateSequence[i+1],gateSequence[i]]=[gateSequence[i],gateSequence[i+1]]; 
      renderGateList();
      container.innerHTML = "<h2>Circuit Diagram</h2>";
      renderCircuit(nQ,gateSequence);
      blochSpheresDiv.innerHTML = "";
      resultsDiv.innerHTML = "";
    } };

    const rm = document.createElement('button');
    rm.textContent = 'Remove';
    rm.className = 'rm';
    rm.onclick = ()=>{ gateSequence.splice(i,1); renderGateList();
      container.innerHTML = "<h2>Circuit Diagram</h2>";
      renderCircuit(nQ,gateSequence);
      blochSpheresDiv.innerHTML = "";
      resultsDiv.innerHTML = "";
    };

    right.appendChild(up);
    right.appendChild(down);
    right.appendChild(rm);

    d.appendChild(left);
    d.appendChild(right);
    gatesListDiv.appendChild(d);
  });
}

// ---------- Quantum ops ----------
function applySingleQubitGate(psi, n, target, U){
  const dim = psi.length;
  const out = Array(dim).fill(0).map(()=>c(0,0));
  for (let i=0;i<dim;i++){
    const bin = i.toString(2).padStart(n, '0');
    const bit = parseInt(bin[target]);
    for (let j=0;j<2;j++){
      const newBin = bin.substring(0,target) + j.toString() + bin.substring(target+1);
      const idx = parseInt(newBin, 2);
      const coeff = U[j][bit];
      out[idx] = math.add(out[idx], math.multiply(coeff, psi[i]));
    }
  }
  return out;
}

function applyCNOT(psi, n, control, target){
  const dim = psi.length;
  const out = Array(dim).fill(0).map(()=>c(0,0));
  for (let i=0;i<dim;i++){
    const bin = i.toString(2).padStart(n, '0');
    if (bin[control] === '1'){
      const flippedBit = bin[target] === '1' ? '0' : '1';
      const newBin = bin.substring(0,target) + flippedBit + bin.substring(target+1);
      const idx = parseInt(newBin, 2);
      out[idx] = math.add(out[idx], psi[i]);
    } else {
      out[i] = math.add(out[i], psi[i]);
    }
  }
  return out;
}

function applyCZ(psi, n, control, target){
  const dim = psi.length;
  const out = Array(dim).fill(0).map(()=>c(0,0));
  for (let i=0;i<dim;i++){
    const bin = i.toString(2).padStart(n, '0');
    const phase = (bin[control]==='1' && bin[target]==='1') ? c(-1,0) : c(1,0);
    out[i] = math.add(out[i], math.multiply(phase, psi[i]));
  }
  return out;
}

function applySWAP(psi, n, a, b){
  if (a===b) return psi.slice();
  const dim = psi.length;
  const out = Array(dim).fill(0).map(()=>c(0,0));
  for (let i=0;i<dim;i++){
    const bin = i.toString(2).padStart(n, '0');
    if (bin[a] === bin[b]){
      out[i] = math.add(out[i], psi[i]);
    } else {
      const swapped = bin.substring(0, Math.min(a,b))
        + (a<b ? bin[b] : bin[a])
        + bin.substring(Math.min(a,b)+1, Math.max(a,b))
        + (a<b ? bin[a] : bin[b])
        + bin.substring(Math.max(a,b)+1);
      const idx = parseInt(swapped, 2);
      out[idx] = math.add(out[idx], psi[i]);
    }
  }
  return out;
}

function applyCCNOT(psi, n, c1, c2, t){
  const dim = psi.length;
  const out = Array(dim).fill(0).map(()=>c(0,0));
  for (let i=0;i<dim;i++){
    const bin = i.toString(2).padStart(n, '0');
    if (bin[c1]==='1' && bin[c2]==='1'){
      const flippedBit = bin[t] === '1' ? '0' : '1';
      const newBin = bin.substring(0,t) + flippedBit + bin.substring(t+1);
      const idx = parseInt(newBin, 2);
      out[idx] = math.add(out[idx], psi[i]);
    } else {
      out[i] = math.add(out[i], psi[i]);
    }
  }
  return out;
}
// ---------- Density matrix & Bloch sphere ----------
function outerProduct(psi){
  const dim = psi.length;
  const rho = Array(dim).fill(0).map(()=>Array(dim).fill(0).map(()=>c(0,0)));
  for (let i=0;i<dim;i++){
    for (let j=0;j<dim;j++){
      rho[i][j] = math.multiply(psi[i], math.conj(psi[j]));
    }
  }
  return rho;
}

function partialTrace(rho, n, target){
  const dim = rho.length;
  let red = [[c(0,0), c(0,0)], [c(0,0), c(0,0)]];
  for (let i=0;i<dim;i++){
    for (let j=0;j<dim;j++){
      const ib = i.toString(2).padStart(n,'0');
      const jb = j.toString(2).padStart(n,'0');
      let equalOther = true;
      for (let k=0;k<n;k++){ if (k!==target && ib[k]!==jb[k]) { equalOther=false; break; } }
      if (equalOther){
        const a = parseInt(ib[target]); const b = parseInt(jb[target]);
        red[a][b] = math.add(red[a][b], rho[i][j]);
      }
    }
  }
  const trace = math.add(red[0][0], red[1][1]);
  return [
    [math.divide(red[0][0], trace), math.divide(red[0][1], trace)],
    [math.divide(red[1][0], trace), math.divide(red[1][1], trace)]
  ];

}

function densityToBloch(red){
  const rho01 = red[0][1];
  const x = 2 * cre(rho01);
  const y = -2 * cim(rho01);
  const z = cre(red[0][0]) - cre(red[1][1]);
  return {x:x, y:y, z:z};
}
function measureQubit(psi, n, target) {
    // psi: state vector array of complex numbers (c(re, im))
    // n: number of qubits
    // target: qubit to measure (0 = leftmost in binary)

    const dim = psi.length;
    let p0 = 0;

    // Compute probability qubit is 0
    for (let i = 0; i < dim; i++) {
        const bin = i.toString(2).padStart(n, '0');
        if (bin[target] === '0') {
            p0 += math.abs(psi[i]) ** 2;
        }
    }

    // Random collapse
    const r = Math.random();
    const outcome = (r < p0) ? 0 : 1;

    // Collapse state
    const newPsi = psi.map((amp, i) => {
        const bin = i.toString(2).padStart(n, '0');
        return (parseInt(bin[target]) === outcome) ? amp : c(0,0);
    });

    // Normalize
    const norm = Math.sqrt(newPsi.reduce((sum, amp) => sum + math.abs(amp)**2, 0));
    for (let i = 0; i < dim; i++) {
        newPsi[i] = math.divide(newPsi[i], norm);
    }

    return { outcome, newPsi };
}

// ---------- Run simulation ----------
function onRun(){
  initState(nQ);
  let psi = stateVec.slice();

  for (const g of gateSequence){
    if (g.type in GATES){
      psi = applySingleQubitGate(psi, nQ, g.params[0], GATES[g.type]);
    }else if(g.type === 'MEASURE'){
      const target = g.params[0];
       // Example: measure qubit 0
      const { outcome, newPsi } = measureQubit(psi, nQ, target);
      psi = newPsi;   // update state vector
      console.log(`Measurement outcome for qubit ${target}:`, outcome);
    }
    else if (g.type === 'Rx'){
      psi = applySingleQubitGate(psi, nQ, g.params[0], Rx(g.angle));
    } else if (g.type === 'Ry'){
      psi = applySingleQubitGate(psi, nQ, g.params[0], Ry(g.angle));
    } else if (g.type === 'Rz'){
      psi = applySingleQubitGate(psi, nQ, g.params[0], Rz(g.angle));
    } else if (g.type === 'Phase'){
      psi = applySingleQubitGate(psi, nQ, g.params[0], Phase(g.angle));
    } else if (g.type === 'CNOT'){
      psi = applyCNOT(psi, nQ, g.params[0], g.params[1]);
    } else if (g.type === 'CZ'){
      psi = applyCZ(psi, nQ, g.params[0], g.params[1]);
    } else if (g.type === 'SWAP'){
      psi = applySWAP(psi, nQ, g.params[0], g.params[1]);
    } else if (g.type === 'CCNOT'){
      psi = applyCCNOT(psi, nQ, g.params[0], g.params[1], g.params[2]);
    }
  }
  stateVec = psi;

  const rho = outerProduct(stateVec);
  const reducedList = [];
  for (let q=0; q<nQ; q++){
    const red = partialTrace(rho, nQ, q);
    reducedList.push(red);
  }

  displayResults(stateVec, rho, reducedList);
  drawAllBloch(reducedList);
}

// ---------- Display & plotting ----------
function displayResults(psi, rho, reducedList){
  resultsDiv.innerHTML = '';
  const dim = psi.length;
  let s = "<div class='result-block'><h3>Final state amplitudes (nonzero)</h3>";
  for (let i = 0; i < dim; i++) {
  const mag = Math.hypot(cre(psi[i]), cim(psi[i]));
  if (mag > 1e-9) {
    const amp = `${cre(psi[i]).toFixed(3)}${cim(psi[i]) >= 0 ? '+' : '-'}${Math.abs(cim(psi[i])).toFixed(3)}j`;
    s += `<div>\\(|${i.toString(2).padStart(nQ,'0')}> : ${amp}\\)</div>`;
  }
}

  s += "</div>";

  s += "<div class='result-block'><h3>Full density matrix ρ</h3>";
  if(rho.length <= 3&& rho[0].length <=3) {
    s += `<div style ="overflow:auto; max-width:100%; max-height: 400px;"><b>$$${formatComplexMatrix(rho)}$$</b><div>`;
  }
  else {
    s += `<div style ="overflow:auto; max-width:100%; max-height: 400px;"><b>${formatMatrixHTML(rho)}</b></div>`;
  }
  s += "</div>";

  for (let q=0;q<reducedList.length;q++){
    s += "<div class='result-block'>";
    s += `<h3>Reduced ρ (qubit ${q})</h3>`;
    const mat = reducedList[q];
    if(mat.length <= 3 && mat[0].length <= 3) {
      s += `<div style ="overflow:auto; max-width:100%; max-height: 400px;">$$${formatComplexMatrix(mat)}$$</div>`;
    }
    else{
      s += `<div style ="overflow:auto; max-width:100%; max-height: 400px;">${formatMatrixHTML(mat)}</div>`;
    }
    s += "</div>";
  }

  resultsDiv.innerHTML = s;
  if(window.MathJax){
    MathJax.typesetPromise();
  }
}

function formatComplexMatrix(mat){
  let latex = "\\begin{bmatrix}\n";
  latex += mat.map(
    row => row.map(
      c => `${cre(c).toFixed(3)}${cim(c) >= 0 ? '+' : ''}${cim(c).toFixed(3)}i`
    ).join(" & ")
  ).join(" \\\\\n");
  latex += "\n\\end{bmatrix}";
  return latex;  
}
function formatMatrixHTML(mat, threshold = 1e-6) {
  return `<table border="1" style="border-collapse: collapse; font-family: monospace;">` +
    mat.map(row => 
      `<tr>` + row.map(c => {
        let val = (Math.hypot(cre(c), cim(c)) < threshold)? "0":`${cre(c).toFixed(2)}${cim(c)>=0?'+':''}${cim(c).toFixed(2)}i`;
         return `<td style = "text-align : center; width : 80px; padding : 2px;">${val}</td>`;
        }).join('') + `</tr>`
    ).join('') +
  `</table>`;
}
function qubitEntropy(x, y, z) {
  // Length of the Bloch vector
  const r = Math.sqrt(x * x + y * y + z * z);

  // Eigenvalues of the density matrix
  const lambda1 = (1 + r) / 2;
  const lambda2 = (1 - r) / 2;

  // Helper for safe log base 2 (avoids log(0) issues)
  function log2Safe(val) {
    return val > 0 ? Math.log(val) / Math.log(2) : 0;
  }

  // Von Neumann entropy
  const S = -(lambda1 * log2Safe(lambda1) + lambda2 * log2Safe(lambda2));
  return S;
}

function drawAllBloch(reducedList) {
  blochSpheresDiv.innerHTML = '';

  for (let q = 0; q < reducedList.length; q++) {
    // wrapper for one qubit (sphere + properties)
    const wrapper = document.createElement('div');
    wrapper.className = 'bloch-wrapper';

    // Bloch sphere div
    const block = document.createElement('div');
    block.id = 'bloch_' + q;
    block.className = 'bloch-canvas';
    wrapper.appendChild(block);
    // Properties div
    const props = document.createElement('div');
    props.className = 'bloch-properties';
    const bloc = densityToBloch(reducedList[q]);
    const x = bloc.x.toFixed(6);
    const y = bloc.y.toFixed(6);
    const z = bloc.z.toFixed(6);
    const e = qubitEntropy(x,y,z).toFixed(3);
    props.innerHTML = `
      <h3>Qubit ${q}</h3>
      <p>Bloch vector:(${bloc.x.toFixed(3)}, ${bloc.y.toFixed(3)}, ${bloc.z.toFixed(3)})</p>
      <p>Entropy(mixedness): ${e}</p>
      <p>Purity: ${(1 + x * x + y * y + z * z) / 2}</p>
      <p>Measurement probabilities(|0>,|1>): ${reducedList[q][0][0]}, ${reducedList[q][1][1]} </p>
    `;
    wrapper.appendChild(props);

    // Add wrapper into output panel
    blochSpheresDiv.appendChild(wrapper);

    // Draw the Bloch sphere
    const r = Math.sqrt(bloc.x**2 + bloc.y**2 + bloc.z**2);
    if(r< 1e-6){
      plotBloch(block.id, {x:0,y:0,z:0},q,true);
    }
    plotBloch(block.id, bloc, q,false);
  }
}

function plotBloch(containerId, bloch, q) {
  const U = 30, V = 30;
  let xs = [], ys = [], zs = [];

  // Sphere coordinates
  for (let i = 0; i <= U; i++) {
    const rowx = [], rowy = [], rowz = [];
    const theta = Math.PI * i / U;
    for (let j = 0; j <= V; j++) {
      const phi = 2 * Math.PI * j / V;
      rowx.push(Math.sin(theta) * Math.cos(phi));
      rowy.push(Math.sin(theta) * Math.sin(phi));
      rowz.push(Math.cos(theta));
    }
    xs.push(rowx); ys.push(rowy); zs.push(rowz);
  }

  // Sphere with wireframe (grid-like mesh)
  const sphere = {
    type: 'surface',
    x: xs, y: ys, z: zs,
    opacity: 0.3,
    colorscale: [[0, 'rgba(228, 246, 253, 0.87)'], [1, 'rgba(248, 200, 244, 1)']],
    showscale: false,
    contours: {
      x: { show: true, color: "#5a56568a", width: 20 },
      y: { show: true, color: "#5a565680", width: 20},
      z: { show: true, color: "#5a565685", width:20 }
    },
    hoverinfo: 'skip'
  };

  // Axes (colored like in your image)
  const axes = [
    { type: 'scatter3d', mode: 'lines', x: [-1, 1], y: [0, 0], z: [0, 0], line: { width: 3, color: 'purple' } ,name:"x-axis"},   // X
    { type: 'scatter3d', mode: 'lines', x: [0, 0], y: [-1, 1], z: [0, 0], line: { width: 3, color: 'purple' } ,name:"y-axis"}, // Y
    { type: 'scatter3d', mode: 'lines', x: [0, 0], y: [0, 0], z: [-1, 1], line: { width: 3, color: 'purple' },name:"z-axis"}     // Z
  ];

  // Qubit vector
  const vx = bloch.x, vy = bloch.y, vz = bloch.z;
  const r = Math.sqrt(vx*vx + vy*vy + vz*vz);
  
  // Arrowhead


  // Basis state labels
  const labels = {
    type: 'scatter3d',
    mode: 'text',
    x: [0, 0, 1.3, -1.3, 0, 0],
    y: [0, 0, 0, 0, 1.3, -1.3],
    z: [1.3, -1.3, 0, 0, 0, 0],
    text: ['z |0⟩', '|1⟩', 'x |+⟩', '|−⟩', 'y |+i⟩', '|−i⟩'],
    textfont: { size: 13, color: '#161618b2' },
    textposition: 'middle center',
    hoverinfo: 'text'
  };

  let traces = [sphere, ...axes, labels];
  if(r<1e-6){
     traces.push({
      type: 'scatter3d',
      mode: 'markers',
      x: [0], y: [0], z: [0],
      marker: { size: 6, color: 'red' },
      name: 'mixed'
    });
  }
  else{
    const stateVector = {
      type: 'scatter3d',
      mode: 'lines+markers',
      x: [0, vx], y: [0, vy], z: [0, vz],
      line: { width: 6, color: '#ff6969ec' },
      marker: { size: 1, color: '#f16464f5' },
      hoverinfo : 'x+y+z',
      name: "qubit"
    };
    const arrowHead = {
      type: 'cone',
      x: [vx], y: [vy], z: [vz],
      u: [vx], v: [vy], w: [vz],   // direction of the vector
      sizemode: 'absolute',
      sizeref: 0.2,
      anchor: 'tip',               // <<< this makes the cone tip sit at (vx,vy,vz)
      colorscale: [[0, 'red'], [1, 'red']],
      showscale: false
    };
    traces.push(stateVector,arrowHead);

  }
  const layout = {
    title: `Qubit ${q}`,
    margin: { l: 0, r: 0, b: 0, t: 30 },
    scene: {
      aspectmode: 'cube',
      xaxis: { range: [-1.3, 1.3], showgrid: false, zeroline: false, showticklabels: false, visible: false },
      yaxis: { range: [-1.3, 1.3], showgrid: false, zeroline: false, showticklabels: false, visible: false },
      zaxis: { range: [-1.3, 1.3], showgrid: false, zeroline: false, showticklabels: false, visible: false },
      camera: { eye: { x: 0.8, y: 0.8, z: 0.8 } }
    }
  };
  Plotly.newPlot(containerId,traces,layout,{displayModeBar : false});
}
document.getElementById("cRun").addEventListener("click", async () => {
  // Build payload from your gates list
  const payload = {
    numQubits: nQ,      // number of qubits from your frontend
    gates: gateSequence // your array of gate objects
  };

  try {
    const res = await fetch("https://qsv-3xax.onrender.com/run", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload)
    });

    const data = await res.json();
    document.getElementById("backendResults").textContent =
      "Backend Counts:\n" + JSON.stringify(data.counts, null, 2) +
      "\n\nQASM:\n" + data.qasm;
  } catch (err) {
    console.error(err);
    document.getElementById("backendResults").textContent = "Error: " + err;
  }
});
 



