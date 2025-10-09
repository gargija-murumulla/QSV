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
const basisInputWrap = document.getElementById('basisInputWrap');
const basisInput = document.getElementById('basisInput');

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

// Custom gate elements
const customMatrixDiv = document.getElementById('customMatrixDiv');
const customCircuitDiv = document.getElementById('customCircuitDiv');
const customControlDiv = document.getElementById('customControlDiv');
const matrixSize = document.getElementById('matrixSize');
const matrixInput = document.getElementById('matrixInput');
const validateMatrix = document.getElementById('validateMatrix');
const matrixValidation = document.getElementById('matrixValidation');
const customGateName = document.getElementById('customGateName');
const circuitGateName = document.getElementById('circuitGateName');
const circuitQubits = document.getElementById('circuitQubits');
const circuitGatesList = document.getElementById('circuitGatesList');
const addCircuitGate = document.getElementById('addCircuitGate');
const validateCircuit = document.getElementById('validateCircuit');
const circuitValidation = document.getElementById('circuitValidation');
const controlGateName = document.getElementById('controlGateName');
const baseGate = document.getElementById('baseGate');
const controlAngleDiv = document.getElementById('controlAngleDiv');
const controlAngle = document.getElementById('controlAngle');
const controlQubitsList = document.getElementById('controlQubitsList');
const controlTarget = document.getElementById('controlTarget');
const controlTarget1 = document.getElementById('controlTarget1');

const validateControl = document.getElementById('validateControl');
const controlValidation = document.getElementById('controlValidation');
const customGatesModal = document.getElementById('customGatesModal');
const closeModal = document.getElementById('closeModal');
const customGatesList = document.getElementById('customGatesList');
const createNewGate = document.getElementById('createNewGate');

const btnAddGate = document.getElementById('btnAddGate');
const btnUndo = document.getElementById('btnUndo');
const btnClearGates = document.getElementById('btnClearGates');
const gatesListDiv = document.getElementById('gatesList');
const btnRun = document.getElementById('btnRun');
const runCircuitBtn = document.getElementById('cRun');
const loadingSpinner = document.getElementById('loading');
const themeToggleBtn = document.getElementById('themeToggle');

const resultsDiv = document.getElementById('results');
const blochSpheresDiv = document.getElementById('blochSpheres');
const container = document.getElementById("quantumCircuit");
  
// ---------- App state ----------
let nQ = 2;
let stateVec = []; // array of math.complex
let gateSequence = []; // each gate object: {type, params, angle?}
let lastData = null;
let customGates = {}; // Store custom gates: {name: {type, matrix, circuit, controls, etc.}}
let currentCustomGate = {
  type: 'CUSTOM_CIRCUIT',
  name: '',
  qubits: nQ,
  gates: []
};
; // Currently being created custom gate
// ---------- Setup handlers ----------
btnSet.addEventListener('click', onSet);
gateType.addEventListener('change', onGateTypeChange);
btnAddGate.addEventListener('click', onAddGate);
btnUndo.addEventListener('click', onUndo);
btnClearGates.addEventListener('click', onClearGates);
btnRun.addEventListener('click', onRun);
// theme toggle
if (themeToggleBtn) {
  initThemeFromStorage();
  themeToggleBtn.addEventListener('click', toggleTheme);
}

// initialize UI
onGateTypeChange();

// Ensure custom gate panels are hidden on page load
customMatrixDiv.classList.add('hidden');
customCircuitDiv.classList.add('hidden');
customControlDiv.classList.add('hidden');
let histogramChart = null;
let targetQubit = null;
let isRunningBackend = false;
// Clear visuals on any input/select change
function clearVisuals(){
  const qDiv = document.getElementById('qsphereDiv');
  if (qDiv){ qDiv.innerHTML = ''; qDiv.classList.add('hidden'); }
  blochSpheresDiv.innerHTML = '';
  resultsDiv.innerHTML = '';
  // clear charts (single-instance policy)
  try {
    if (window.stateChart) {
      window.stateChart.data.labels = [];
      window.stateChart.data.datasets[0].data = [];
      window.stateChart.update();
    }
  } catch(e) {}
  try {
    if (typeof histogramChart !== 'undefined' && histogramChart) {
      histogramChart.data.labels = [];
      histogramChart.data.datasets[0].data = [];
      histogramChart.update();
    }
  } catch(e) {}
}
document.addEventListener('change', (e)=>{
  const t = e.target;
  if (t && (t.tagName === 'INPUT' || t.tagName === 'SELECT')) {
    clearVisuals();
  }
});
document.addEventListener('input', (e)=>{
  const t = e.target;
  if (t && (t.tagName === 'INPUT' || t.tagName === 'SELECT')) {
    clearVisuals();
  }
});
// ---------- Functions ----------
blochSpheresDiv.innerHTML = "<div class = grid > <b><h2 style= text-align:'centre'>Qubit</h2></b> <p>The fundemental unit of quantum information, serving as the quantum equvivalent of a classical computer's bit. A qubit can have states 0, 1, 0/1(superposition). </p></div>"
function drawHistogram(counts) {
    const labels = Object.keys(counts);
    const values = Object.values(counts);
    const canvas = document.getElementById('histogram');
    const ctx = canvas.getContext('2d', { willReadFrequently: true });

    // Reuse existing chart if present
    const existing = Chart.getChart(canvas);
    if (existing) {
      histogramChart = existing;
      histogramChart.data.labels = labels;
      histogramChart.data.datasets[0].data = values;
      histogramChart.update();
      return;
    }

    histogramChart = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: labels,
            datasets: [{
                label: 'Measurement Counts',
                data: values,
                backgroundColor: 'rgba(149, 206, 243, 0.91)',
                borderColor: 'rgb(7, 74, 119)',
                borderWidth: 1
            }]
        },
        options: {
            responsive: true,
            plugins: {
                legend: {
                    labels: { color: '#000000' }
                }
            },
            scales: {
                x: {
                    title: { display: true, text: 'Bitstring Outcome', color: '#000000' },
                    ticks: { color: '#000000' },
                    grid: { color: 'rgba(0,0,0,0.1)' }
                },
                y: {
                    beginAtZero: true,
                    title: { display: true, text: 'Counts', color: '#000000' },
                    ticks: { color: '#000000' },
                    grid: { color: 'rgba(0,0,0,0.1)' }
                }
            }
        }
    });
}

let firstQubit = false;
function onSet(){
  nQ = parseInt(numQInput.value);
  // Switch input UI: dropdown for <=10, text input for >10
  if (nQ > 10) {
    basisSelect.classList.add('hidden');
    basisInputWrap.classList.remove('hidden');
  } else {
    basisInputWrap.classList.add('hidden');
    basisSelect.classList.remove('hidden');
  }
  if ((nQ >=1 && nQ <=5)) { 
    populateBasis(nQ);
    populateQubitSelectors(nQ);
    const initIndex = 0;
    initStates = getInitStates(initIndex, nQ);
    console.log("👉 Default initial state:", initStates);
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
    customQubitContainer.classList.add("hidden");
    btnRun.classList.remove("hidden");
    stateBarChart.classList.add("hidden");
  }
  else{
    populateBasis(nQ);
    populateQubitSelectors(nQ);
    const initIndex = 0;
    initStates = getInitStates(initIndex, nQ);
    console.log("👉 Default initial state:", initStates);
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

    customQubitContainer.classList.remove("hidden");
    btnRun.classList.add("hidden");
    const qubitSelect = document.getElementById("qubitSelect");
    qubitSelect.innerHTML = "";
    for (let i = 0; i < nQ; i++) {
      const opt = document.createElement("option");
      opt.value = i;
      opt.textContent = "q" + i;
      if(i===0){
        opt.selected = true;
      }
      qubitSelect.appendChild(opt);
    }
    

  }
}
function populateBasis(n){
  if (n > 10) return; // skip generating huge dropdown
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
  const sels = [targetQ, controlQ, controlTarget,controlTarget1, targetQ2, swapA, swapB, cc_c1, cc_c2, cc_t];
  sels.forEach(s => s.innerHTML = '');
  for (let i=0;i<n;i++){
    const opt = (id)=>{ const o=document.createElement('option'); o.value=i; o.text='q'+i; return o; };
    sels.forEach(s => s.appendChild(opt()));
  }
}
let initStates =[];
function initState(nQ){
  const dim = 1<<nQ;
  let initIndex = 0;
  if (nQ > 10) {
    const raw = (basisInput.value || '').trim();
    const valid = raw && raw.length === nQ && /^[01]+$/.test(raw);
    if (!valid) {
      alert(`Please enter a ${nQ}-bit binary string (0/1).`);
      initIndex = 0;
      basisInput.value = '0'.repeat(nQ);
    } else {
      initIndex = parseInt(raw, 2);
    }
  } else {
    initIndex = parseInt(basisSelect.value || "0", 2);
  }
  stateVec = Array(dim).fill(0).map(()=>c(0,0));
  stateVec[initIndex] = c(1,0);
  initStates = getInitStates(initIndex,nQ);
}

function plotQSphere(divId, stateVec) {
  const colors = getPlotColors();
  const nQ = Math.log2(stateVec.length);
  const spikeTraces = [];
  const tipTraces = [];
  const latitudeTraces = [];
  const labelX = [];
  const labelY = [];
  const labelZ = [];
  const labelText = [];

  const coords = [];
  const phases = [];
  const probs = [];
  const arcTraces = [];

  // --- Compute spike positions ---
  for (let i = 0; i < stateVec.length; i++) {
    const amp = stateVec[i];
    const re = amp.re, im = amp.im;
    const prob = re*re + im*im;
    const phase = Math.atan2(im, re);
    const weightStr = i.toString(2).padStart(nQ,'0');

    // --- evenly distribute states on sphere ---
    const hamming = weightStr.split('').filter(q => q==='1').length;
    const theta = (hamming / nQ) * Math.PI;           // latitude by Hamming weight
    const phi = 2 * Math.PI * i / stateVec.length;   // evenly around longitude

    const r = 1.0;
    const x = r * Math.sin(theta) * Math.cos(phi);
    const y = r * Math.sin(theta) * Math.sin(phi);
    const z = r * Math.cos(theta);

    coords.push([x, y, z]);
    phases.push(phase);
    probs.push(prob);

    labelX.push(x);
    labelY.push(y);
    labelZ.push(z);
    labelText.push(`|${weightStr}⟩`);

    // radial line
    spikeTraces.push({
      type:"scatter3d",
      mode:"lines",
      x:[0, x], y:[0, y], z:[0, z],
      line:{color:`hsl(${(phase*180/Math.PI+360)%360},80%,50%)`, width:1 + 8*prob},
      opacity:0.8,
      hoverinfo:"skip",
      showlegend : false
    });

    // tip marker for hover text
    tipTraces.push({
      type: "scatter3d",
      mode: "markers",
      x: [x], y: [y], z: [z],
      marker: {size: 5 + 20*prob, color:`hsl(${(phase*180/Math.PI+360)%360},80%,40%)`},
      text: `|${weightStr}⟩<br>amp=${re.toFixed(2)} + ${im.toFixed(2)}i<br>P=${prob.toFixed(2)}<br>phase=${phase.toFixed(2)}`,
      hoverinfo: "text",
      showlegend : false
    });
  }
  // --- Latitude circles ---
  const SPHERE_POINTS = 60;
  for (let k = 0; k <= nQ; k++) {
    const theta = (k / nQ) * Math.PI;
    const latX = [], latY = [], latZ = [];
    for (let p = 0; p <= SPHERE_POINTS; p++) {
      const phi = (p / SPHERE_POINTS) * 2 * Math.PI;
      latX.push(Math.sin(theta)*Math.cos(phi));
      latY.push(Math.sin(theta)*Math.sin(phi));
      latZ.push(Math.cos(theta));
    }
    latitudeTraces.push({
      type:"scatter3d",
      mode:"lines",
      x:latX, y:latY, z:latZ,
      line:{color:"gray", width:1},
      opacity:0.2,
      hoverinfo:"skip",
      showlegend:false
    });
  }

  // --- Transparent sphere ---
  const U = 30, V = 30;
  const xs = [], ys = [], zs = [];
  for (let i = 0; i <= U; i++) {
    const theta = Math.PI * i / U;
    const rowX = [], rowY = [], rowZ = [];
    for (let j = 0; j <= V; j++) {
      const phi = 2*Math.PI*j/V;
      rowX.push(Math.sin(theta)*Math.cos(phi));
      rowY.push(Math.sin(theta)*Math.sin(phi));
      rowZ.push(Math.cos(theta));
    }
    xs.push(rowX); ys.push(rowY); zs.push(rowZ);
  }

  const sphereSurface = {
    type:'surface', x:xs, y:ys, z:zs,
    opacity:0.2,
    colorscale:[[0,'rgba(228,246,253,0.87)'], [1,'rgba(248,200,244,0.5)']],
    showscale:false,
    contours: {
      x: { show: true, color: "#5a56568a", width: 20 },
      y: { show: true, color: "#5a565680", width: 20},
      z: { show: true, color: "#5a565685", width:20 }
    },
    hoverinfo:'skip',
    showlegend:false
  };

  const labelTraces = {
    type:"scatter3d",
    mode:"text",
    x:labelX, y:labelY, z:labelZ,
    text:labelText,
    textposition:"top center",
    textfont:{size:12, color:"#333333"},
    hoverinfo:"skip",
    showlegend:false
  };
  const layout = {
    title:"Q-Sphere <br> size(Dot)-> probability<br> color(Dot)->Phase",
    margin:{l:0,r:0,b:0,t:30},
    paper_bgcolor: colors.paper,
    plot_bgcolor: colors.plot,
    scene:{
      bgcolor: colors.scene,
      aspectmode:'cube',
      xaxis:{range:[-1.3,1.3],showgrid:false,zeroline:false,showticklabels:false,visible:false},
      yaxis:{range:[-1.3,1.3],showgrid:false,zeroline:false,showticklabels:false,visible:false},
      zaxis:{range:[-1.3,1.3],showgrid:false,zeroline:false,showticklabels:false,visible:false},
      camera:{eye:{x:0.8,y:0.8,z:0.8}}
    },
  };

  Plotly.newPlot(divId, [sphereSurface, ...latitudeTraces, ...spikeTraces, ...tipTraces, labelTraces], layout);
}

function getInitStates(initIndex,nQ) {
  const initStates =[];
  for (let i = nQ - 1; i >= 0; i--) {
    initStates.push((initIndex >> i) & 1);
  }
  return initStates;
}


 function onGateTypeChange(){
  const type = gateType.value;
  // hide all
  singleTargetDiv.classList.add('hidden');
  cnotDiv.classList.add('hidden');
  swapDiv.classList.add('hidden');
  ccnotDiv.classList.add('hidden');
  angleDiv.classList.add('hidden');
  customMatrixDiv.classList.add('hidden');
  customCircuitDiv.classList.add('hidden');
  customControlDiv.classList.add('hidden');

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
  if (type === 'CUSTOM_MATRIX'){
    customMatrixDiv.classList.remove('hidden');
    setupMatrixInput();
  }
  if (type === 'CUSTOM_CIRCUIT'){
    customCircuitDiv.classList.remove('hidden');
    setupCircuitBuilder();
  }
  if (type === 'CUSTOM_CONTROL'){
    customControlDiv.classList.remove('hidden');
    setupControlBuilder();
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
  } else if (type === 'CUSTOM_MATRIX'){
    if (!validateCustomMatrix()) return;
    const name = customGateName.value.trim();
    const target = parseInt(controlTarget.value);
    gate = customGates[name];
  } else if (type === 'CUSTOM_CIRCUIT'){
    if (!validateCustomCircuit()) return;
    const name = circuitGateName.value.trim();
    gate = customGates[name];
  } else if (type === 'CUSTOM_CONTROL'){
    if (!validateCustomControl()) return;
    const name = controlGateName.value.trim();
    const selectedControls = Array.from(controlQubitsList.querySelectorAll('input[type="checkbox"]:checked'))
      .map(cb => parseInt(cb.value));
    const target = parseInt(controlTarget.value);
    gate = customGates[name];
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

    // Custom gates
    if (g.type === "CUSTOM") {
      const customGate = customGates[g.name];
      if (customGate) {
        const rect = document.createElementNS(svgNS, "rect");
        rect.setAttribute("x", x - 30);
        rect.setAttribute("y", 30 - 25);
        rect.setAttribute("width", 60);
        rect.setAttribute("height", numQubits * qheight + 50);
        rect.setAttribute("fill", "#e8f4fd");
        rect.setAttribute("stroke", "#2196F3");
        rect.setAttribute("stroke-width", "2");
        svg.appendChild(rect);

        const label = document.createElementNS(svgNS, "text");
        label.setAttribute("x", x);
        label.setAttribute("y", 30 + (numQubits * qheight) / 2);
        label.setAttribute("text-anchor", "middle");
        label.setAttribute("dominant-baseline", "middle");
        label.setAttribute("font-size", "12");
        label.setAttribute("font-weight", "bold");
        label.textContent = g.name;
        svg.appendChild(label);

        // Add type indicator
        const typeLabel = document.createElementNS(svgNS, "text");
        typeLabel.setAttribute("x", x);
        typeLabel.setAttribute("y", 30 + (numQubits * qheight) / 2 + 15);
        typeLabel.setAttribute("text-anchor", "middle");
        typeLabel.setAttribute("dominant-baseline", "middle");
        typeLabel.setAttribute("font-size", "10");
        typeLabel.setAttribute("fill", "#666");
        typeLabel.textContent = customGate.type.replace('CUSTOM_', '');
        svg.appendChild(typeLabel);
      }
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
      } else if (g.type === "SWAP"){
        isTarget = (q === g.params[0] || q === g.params[1]);
      } else if (g.type === "CUSTOM") {
        // Custom gates affect all qubits for now
        isTarget = true;
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
    let desc = `${i+1}.`;
    if (g.type === 'CUSTOM') {
      
      desc += ` ${g.name}`;
      const customGate = customGates[g.name];
      if (customGate) {
        desc = `${i+1}. ${g.name} (Custom ${customGate.type.replace('CUSTOM_', '')})`;
        if (g.params?.length) {
          desc += ` (${g.params.join(',')})`;
        }
      } else {
        desc = `${i+1}. ${g.name} (Custom Gate - Not Found)`;
      }
    } else {
      desc += ` ${g.type}`;
      if (g.params?.length){
        desc += ` (${g.params.join(',')})`;
      }
      if (g.angle !== undefined){
        const deg = (g.angle*180/Math.PI).toFixed(2);
        desc += `, ${deg}°`;
      }
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

// Custom gate application functions
function applyCustomGate(psi, n, gate) {
  const customGate = customGates[gate.name];
  if (!customGate) {
    console.error(`Custom gate "${gate.name}" not found`);
    return psi;
  }
  
  console.log(`Applying custom gate: ${gate.name} (${customGate.type})`);
  console.log('Input state vector:', psi.map((amp, i) => `${i}: ${cre(amp).toFixed(3)}+${cim(amp).toFixed(3)}i`));
  
  let result;
  switch (customGate.customType) {
    case 'CUSTOM_MATRIX':
      result = applyCustomMatrixGate(psi, n, customGate, gate.params[0]);
      break;
    case 'CUSTOM_CIRCUIT':
      result = applyCustomCircuitGate(psi, n, customGate);
      break;
    case 'CUSTOM_CONTROL':
      result = applyCustomControlGate(psi, n, customGate);
      break;
    default:
      console.error(`Unknown custom gate type: ${customGate.type}`);
      return psi;
  }
  
  console.log('Output state vector:', result.map((amp, i) => `${i}: ${cre(amp).toFixed(3)}+${cim(amp).toFixed(3)}i`));
  
  // Validate the result
  for (let i = 0; i < result.length; i++) {
    if (isNaN(cre(result[i])) || isNaN(cim(result[i]))) {
      console.error(`NaN in result at state ${i}:`, result[i]);
      return psi; // Return original state if NaN detected
    }
  }
  
  return result;
}

function applyCustomMatrixGate(psi, n, customGate, target) {
  const matrix = customGate.matrix;
  const size = customGate.size;
  const numQubits = Math.log2(size);
  
  console.log(`Applying custom matrix gate: ${numQubits} qubits, target: ${target}`);
  console.log('Matrix:', matrix);
  
  if (numQubits === 1) {
    // Single qubit gate - use existing function
    return applySingleQubitGate(psi, n, target, matrix);
  } else {
    // Multi-qubit gate - simplified implementation
    // For now, just apply as single qubit to the target
    console.warn(`Multi-qubit custom matrix gate not fully implemented. Applying as single qubit to target ${target}`);
    return applySingleQubitGate(psi, n, target, matrix);
  }
}
function applyCustomCircuitGate(psi, n, customGate) {
  let result = psi.slice();
  
  // Apply each gate in the circuit
  for (const gate of customGate.subGates) {
    if (gate.type in GATES) {
      result = applySingleQubitGate(result, n, gate.params[0], GATES[gate.type]);
    } else if (gate.type === 'Rx') {
      result = applySingleQubitGate(result, n, gate.params[0], Rx(gate.angle));
    } else if (gate.type === 'Ry') {
      result = applySingleQubitGate(result, n, gate.params[0], Ry(gate.angle));
    } else if (gate.type === 'Rz') {
      result = applySingleQubitGate(result, n, gate.params[0], Rz(gate.angle));
    } else if (gate.type === 'Phase') {
      result = applySingleQubitGate(result, n, gate.params[0], Phase(gate.angle));
    } else if (gate.type === 'CNOT') {
      result = applyCNOT(result, n, gate.params[1], gate.params[0]);
    } else if (gate.type === 'CZ') {
      result = applyCZ(result, n, gate.params[1], gate.params[0]);
    } else if (gate.type === 'SWAP') {
      result = applySWAP(result, n, gate.params[1], gate.params[0]);
    }else if (gate.type === 'CCNOT') {
      result = applyCCNOT(result, n, gate.params[1], gate.params[0]);
    }
  }
  
  return result;
}
function applyCustomControlGate(psi, n, customGate) {
  // Extract information from new structure
  if (customGate.customType !== 'CUSTOM_CONTROL') return psi;

  const control = customGate.params[0];
  const target = customGate.params[1];
  const subGate = customGate.subGates[0];
  const baseGate = subGate.type;
  const angle = subGate.angle || customGate.angle || 0;

  const result = psi.map(() => c(0,0));

  // Choose gate matrix
  let gateMatrix;
  if (baseGate === 'X') gateMatrix = [[c(0,0), c(1,0)], [c(1,0), c(0,0)]];
  else if (baseGate === 'Y') gateMatrix = [[c(0,0), c(0,-1)], [c(0,1), c(0,0)]];
  else if (baseGate === 'Z') gateMatrix = [[c(1,0), c(0,0)], [c(0,0), c(-1,0)]];
  else if (baseGate === 'H') {
    const s2 = 1 / Math.sqrt(2);
    gateMatrix = [[c(s2,0), c(s2,0)], [c(s2,0), c(-s2,0)]];
  } else if (baseGate === 'Rx') gateMatrix = Rx(angle);
  else if (baseGate === 'Ry') gateMatrix = Ry(angle);
  else if (baseGate === 'Rz') gateMatrix = Rz(angle);
  else if (baseGate === 'Phase') gateMatrix = Phase(angle);
  else return psi; // Unknown subgate type

  // Apply controlled gate logic
  for (let i = 0; i < psi.length; i++) {
    const bin = i.toString(2).padStart(n, '0');
    const controlBit = bin[control];
    const targetBit = bin[target];

    // If control not active → just copy
    if (controlBit !== '1') {
      result[i] = psi[i];
      continue;
    }

    // Flip target bit to get partner state index
    const flippedBin =
      bin.substring(0, target) +
      (targetBit === '0' ? '1' : '0') +
      bin.substring(target + 1);
    const j = parseInt(flippedBin, 2);

    const a0 = psi[i];
    const a1 = psi[j];

    // Apply gate matrix to amplitudes
    const new0 = math.add(
      math.multiply(gateMatrix[0][0], a0),
      math.multiply(gateMatrix[0][1], a1)
    );
    const new1 = math.add(
      math.multiply(gateMatrix[1][0], a0),
      math.multiply(gateMatrix[1][1], a1)
    );

    result[i] = new0;
    result[j] = new1;
  }

  return result;
}


// function applyCustomControlGate(psi, n, customGate) {
//   const { baseGate, controls, target, angle } = customGate;
//   const result = Array(psi.length).fill(null).map(() => c(0,0));

//   // Choose the gate matrix
//   let gateMatrix;
//   if (baseGate === 'X') gateMatrix = [[c(0,0), c(1,0)], [c(1,0), c(0,0)]];
//   else if (baseGate === 'Y') gateMatrix = [[c(0,0), c(0,-1)], [c(0,1), c(0,0)]];
//   else if (baseGate === 'Z') gateMatrix = [[c(1,0), c(0,0)], [c(0,0), c(-1,0)]];
//   else if (baseGate === 'H') {
//     const s2 = 1/Math.sqrt(2);
//     gateMatrix = [[c(s2,0), c(s2,0)], [c(s2,0), c(-s2,0)]];
//   } 
//   else if (baseGate === 'Rx') gateMatrix = Rx(angle || 0);
//   else if (baseGate === 'Ry') gateMatrix = Ry(angle || 0);
//   else if (baseGate === 'Rz') gateMatrix = Rz(angle || 0);
//   else if (baseGate === 'Phase') gateMatrix = Phase(angle || 0);
//   else return psi; // unknown gate

//   // Loop over all basis states
//   for (let i = 0; i < psi.length; i++) {
//     const bin = i.toString(2).padStart(n, '0');

//     // Check if controls are all 1
//     const allControlsOne = controls.every(c => bin[c] === '1');
//     if (!allControlsOne) {
//       // If controls not satisfied → copy amplitude
//       result[i] = math.add(result[i], psi[i]);
//       continue;
//     }

//     // Find partner index: same as i except target flipped
//     const flippedBin = bin.substring(0,target) + (bin[target] === '0' ? '1' : '0') + bin.substring(target+1);
//     const j = parseInt(flippedBin, 2);

//     // Only handle each pair once (when target=0)
//     if (bin[target] === '0') {
//       const a0 = psi[i] || c(0,0);
//       const a1 = psi[j] || c(0,0);

//       // Apply 2x2 gate matrix
//       const new0 = math.add(math.multiply(gateMatrix[0][0], a0), math.multiply(gateMatrix[0][1], a1));
//       const new1 = math.add(math.multiply(gateMatrix[1][0], a0), math.multiply(gateMatrix[1][1], a1));

//       result[i] = math.add(result[i], new0);
//       result[j] = math.add(result[j], new1);
//     }
//   }

//   return result;
// }

function applyMultiQubitGate(psi, n, matrix, targets) {
  const k = targets.length;  // number of target qubits
  const dim = psi.length;    // = 2^n
  const subDim = 1 << k;     // size of the gate (2^k)
  
  if (matrix.length !== subDim) {
    throw new Error(`Matrix size ${matrix.length} does not match target qubits (${k})`);
  }

  const result = Array(dim).fill(null).map(() => c(0,0));
  const targetMask = targets.reduce((acc, t) => acc | (1 << (n-1-t)), 0); // bitmask for target qubits

  // Loop through all basis states
  for (let base = 0; base < dim; base++) {
    // Extract the index of just the target qubits
    let targetIndex = 0;
    for (let j = 0; j < k; j++) {
      const bitPos = n - 1 - targets[j];
      const bit = (base >> bitPos) & 1;
      targetIndex |= (bit << (k - 1 - j));
    }

    // Build the input subvector
    const subvec = Array(subDim).fill(c(0,0));
    subvec[targetIndex] = psi[base];

    // Apply the matrix
    for (let row = 0; row < subDim; row++) {
      let accum = c(0,0);
      for (let col = 0; col < subDim; col++) {
        accum = math.add(accum, math.multiply(matrix[row][col], subvec[col]));
      }

      // Write back to correct global index
      let newIndex = base;
      for (let j = 0; j < k; j++) {
        const bit = (row >> (k - 1 - j)) & 1;
        const bitPos = n - 1 - targets[j];
        if (bit) {
          newIndex |= (1 << bitPos);
        } else {
          newIndex &= ~(1 << bitPos);
        }
      }

      result[newIndex] = math.add(result[newIndex], accum);
    }
  }

  return result;
}

// function applyMultiQubitGate(psi, n, matrix, targets) {
//   // This is a simplified implementation for multi-qubit gates
//   // In a full implementation, you'd need to handle the tensor product properly
//   const dim = psi.length;
//   const out = Array(dim).fill(0).map(() => c(0, 0));
  
//   // For now, just apply as single qubit to first target
//   if (targets.length > 0) {
//     return applySingleQubitGate(psi, n, targets[0], matrix);
//   }
  
//   return out;
// }

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
    if(nQ<=5){
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
      } else if (g.type === 'CUSTOM'){
        psi = applyCustomGate(psi, nQ, g);
        
        // Debug: Check for NaN values after custom gate application
        for (let i = 0; i < psi.length; i++) {
          if (isNaN(cre(psi[i])) || isNaN(cim(psi[i]))) {
            console.error(`NaN detected after applying custom gate ${g.name} at state ${i}:`, psi[i]);
            console.error('Gate details:', g);
            console.error('Custom gate definition:', customGates[g.name]);
          }
        }
      }
    }
    qsphereDiv.classList.remove("hidden");
    stateVec = psi;
    
    // Debug: Log final state vector
    console.log('Final state vector after all gates:');
    stateVec.forEach((amp, i) => {
      const prob = cre(amp) * cre(amp) + cim(amp) * cim(amp);
      if (prob > 1e-6) {
        console.log(`|${i.toString(2).padStart(nQ, '0')}⟩: ${cre(amp).toFixed(3)}+${cim(amp).toFixed(3)}i (P=${prob.toFixed(3)})`);
      }
    });
    
    plotQSphere("qsphereDiv", stateVec);
    const rho = outerProduct(stateVec);
    const reducedList = [];
    for (let q=0; q<nQ; q++){
      const red = partialTrace(rho, nQ, q);
      reducedList.push(red);
    }
    console.log(stateVec);
    console.log(rho);
    console.log(reducedList);
    displayResults(stateVec, rho, reducedList);
    drawAllBloch(reducedList);
  }
  else{
    alert("for more than 5 qubits we use tomography ");
  }
}

function drawStateBarGraph(stateVec) {
  // Compute probabilities (magnitude squared)
  const probs = stateVec.map(a => a.re * a.re + a.im * a.im);

  // Labels as binary states
  const labels = probs.map((_, i) => i.toString(2).padStart(Math.log2(probs.length), "0"));

  if (window.stateChart) {
    window.stateChart.data.labels = labels;
    window.stateChart.data.datasets[0].data = probs;

    // call Chart.js update method
    window.stateChart.update()
    return;
  }

  const ctx = document.getElementById("stateBarChart").getContext("2d");
  window.stateChart = new Chart(ctx, {
    type: "bar",
    data: {
      labels: labels,
      datasets: [{
        label: "Probability",
        data: probs,
        backgroundColor: "rgba(54, 162, 235, 0.7)",
        borderColor: "rgba(54, 162, 235, 1)",
        borderWidth: 1
      }]
    },
    options: {
      responsive: true,
      plugins: {
        legend: {
          labels: { color: '#000000' }
        }
      },
      scales: {
        x: {
          ticks: { color: '#000000' },
          grid: { color: 'rgba(0,0,0,0.1)' }
        },
        y: {
          beginAtZero: true,
          max: 1,
          ticks: { color: '#000000' },
          grid: { color: 'rgba(0,0,0,0.1)' }
        }
      }
    }
  });
}
// ---------- Display & plotting ----------
function displayResults(psi, rho, reducedList){
  drawStateBarGraph(stateVec);

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
function toComplex(c) {
  if (c && typeof c === "object" && "real" in c && "imag" in c) {
    return math.complex(c.real, c.imag);
  }
  return c; // already a math.js complex or number
}

function formatComplexMatrix(mat) {
  let latex = "\\begin{bmatrix}\n";
  latex += mat.map(
    row => row.map(c => {
      const cc = toComplex(c);
      return `${cc.re.toFixed(3)}${cc.im >= 0 ? '+' : ''}${cc.im.toFixed(3)}i`;
    }).join(" & ")
  ).join(" \\\\\n");
  latex += "\n\\end{bmatrix}";
  return latex;  
}
function formatMatrixHTML(mat, threshold = 1e-6) {
  if (!mat || !Array.isArray(mat)) {
    return "<i>Invalid matrix</i>";
  }

  return `<table border="1" style="border-collapse: collapse; font-family: monospace;">` +
    mat.map(row => 
      `<tr>` + row.map(c => {
        const cc = toComplex(c); // ✅ ensure proper math.js complex
        let val = (Math.hypot(cc.re, cc.im) < threshold)
          ? "0"
          : `${cc.re.toFixed(2)}${cc.im >= 0 ? '+' : ''}${cc.im.toFixed(2)}i`;
        return `<td style="text-align:center; width:80px; padding:2px;">${val}</td>`;
      }).join('') + `</tr>`
    ).join('') +
  `</table>`;
}
function formatComplex(entry, digits = 3) {
  const re = Number(entry.re).toFixed(digits);
  const im = Number(entry.im).toFixed(digits);
  return `${re}${im >= 0 ? " + " : " - "}${Math.abs(im)}i`;
}

function formatReal(entry, digits = 3) {
  return Number(entry.re).toFixed(digits);
}

function qubitEntropy(x, y, z) {
  // Length of the Bloch vector
  const r = Math.sqrt(x * x + y * y + z * z);

  // Eigenvalues of the density matrix
  const lambda1 = (1 + r) / 2;
  const lambda2 = (1 - r) / 2;

  // Helper for safe log base 2 (avoids log(0) issues)

  // Von Neumann entropy
  const S = -(lambda1 * log2Safe(lambda1) + lambda2 * log2Safe(lambda2));
  return S;
}
function log2Safe(val) {
  return val > 0 ? Math.log(val) / Math.log(2) : 0;
}
function roundVal(val, digits=3) {
  const num = Number(val.toFixed(digits));
  return Math.abs(num) < 1e-12 ? 0 : num;
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
    const e =cleanFixed(qubitEntropy(x,y,z).toFixed(3));
    props.innerHTML = `
      <h3>Qubit ${q}</h3>
      <p>Bloch vector:(${bloc.x.toFixed(3)}, ${bloc.y.toFixed(3)}, ${bloc.z.toFixed(3)})</p>
      <p>Entropy(mixedness): ${e}</p>
      <p>Purity: ${((1 + x * x + y * y + z * z) / 2).toFixed(3)}</p>
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
  const colors = getPlotColors();
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
    paper_bgcolor: colors.paper,
    plot_bgcolor: colors.plot,
    scene: {
      bgcolor: colors.scene,
      aspectmode: 'cube',
      xaxis: { range: [-1.3, 1.3], showgrid: false, zeroline: false, showticklabels: false, visible: false },
      yaxis: { range: [-1.3, 1.3], showgrid: false, zeroline: false, showticklabels: false, visible: false },
      zaxis: { range: [-1.3, 1.3], showgrid: false, zeroline: false, showticklabels: false, visible: false },
      camera: { eye: { x: 0.8, y: 0.8, z: 0.8 } }
    }
  };
  Plotly.newPlot(containerId,traces,layout,{displayModeBar : false});
}
function cleanFixed(val, digits = 3) {
  const num = Number(val);              // force conversion to number
  const rounded = Number(num.toFixed(digits));
  return rounded === 0 ? 0 : rounded;   // remove negative zero
}


document.getElementById("qubitSelect").addEventListener("change", (e) => {
  const q = parseInt(e.target.value);
  showQubitInfo(q);
});
basisSelect.addEventListener("change", () => {
  const initIndex = parseInt(basisSelect.value, 2);
  initStates = getInitStates(initIndex, nQ);
  console.log("👉 Updated initial states:", initStates);
});
if (basisInput) {
  basisInput.addEventListener('input', ()=>{
    const raw = basisInput.value.trim();
    if (raw.length > nQ) basisInput.value = raw.slice(0, nQ);
  });
}
document.getElementById("cRun").addEventListener("click", async () => {
  // legacy direct handler retained but wrapped above
  // Build payload from your gates list
  
  try {
    
    if(nQ<6){
      const payload = {
        numQubits: nQ,      // number of qubits from your frontend
        gates: gateSequence, // your array of gate objects
        initialStates : initStates, // initial states if any
      };

      console.log("👉 Sending to backend:", payload);
      const res = await fetch("https://qsv-3xax.onrender.com/run", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload)
      });
      const data = await res.json();
      document.getElementById("backendResults").textContent =
        "Backend Counts:\n" + JSON.stringify(data.counts, null, 2) +
        "\n\nQASM:\n" + data.qasm;
      drawHistogram(data.counts);
    }
    else{
      alert("we used single qubit Quantum State tomography for qubits more than 5. The results may have some eroor due to little noise ");
      targetQubit = parseInt(document.getElementById("qubitSelect").value);
      if (isNaN(targetQubit)) {
        alert("Please select a valid target qubit!");
        return;
      }
      console.log(targetQubit);
      const payload = {
        numQubits: nQ,      // number of qubits from your frontend
        gates: gateSequence, // your array of gate objects
        initialStates : initStates, // initial states if any
        targetQubit: targetQubit 
      };
  
      console.log("👉 Sending to backend:", payload);
      const res = await fetch("https://qsv-3xax.onrender.com/run", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload)
      });
      const data = await res.json();
      
      console.log(document.getElementById("qubitSelect").innerHTML);
      console.log("👉 Backend returned:", data);
      // show first qubit immediately

      let x = roundVal(data.blochs.x,3);
      let y = roundVal(data.blochs.y,3);
      let z = roundVal(data.blochs.z,3);
      document.getElementById("backendResults").textContent =
        `Bloch : (${x}, ${y}, ${z})`+
        "\n\nRhos:\n" + formatComplexMatrix(data.rho);
      MathJax.typeset();
      // clear any previous bloch wrappers for tomography case
      blochSpheresDiv.innerHTML = '';
      const wrapper = document.createElement('div');
      wrapper.className = 'bloch-wrapper';

      // Bloch sphere div
      const block = document.createElement('div');
      block.id = 'bloch_' + targetQubit;
      block.className = 'bloch-canvas';
      wrapper.appendChild(block);
      // Properties div
      const props = document.createElement('div');
      props.className = 'bloch-properties';
      const e = cleanFixed(qubitEntropy(x,y,z).toFixed(3));
      props.innerHTML = `
        <h3>Qubit ${targetQubit}</h3>
        <p>Bloch vector:(${x}, ${y}, ${z})</p>
        <p>Entropy(mixedness): ${e}</p>
        <p>Purity: ${(1 + x * x + y * y + z * z) / 2}</p>
        <p>Measurement probabilities(|0>,|1>): 
           ${formatReal(data.rho[0][0], 3)}, 
           ${formatReal(data.rho[1][1], 3)}
        </p>
      `;
      wrapper.appendChild(props);

      // Add wrapper into output panel
      blochSpheresDiv.appendChild(wrapper);
      x = data.blochs.x.toFixed(3);
      y = data.blochs.y.toFixed(3);
      z = data.blochs.z.toFixed(3);
        // Draw the Bloch sphere
      const r = Math.sqrt(x**2 + y**2 + z**2);
      if(r< 1e-6){
        plotBloch(block.id, {x:0,y:0,z:0}, targetQubit,true);
      }
      plotBloch(block.id,{x,y,z},targetQubit,true);
    }

  }
  catch (err) {
      console.error(err);
      document.getElementById("backendResults").textContent = "Error: " + err;
  }
});
// ---------- Custom Gate Functions ----------

// Matrix-based custom gate setup
function setupMatrixInput() {
  const size = parseInt(matrixSize.value);
  matrixInput.innerHTML = '';
  
  const container = document.createElement('div');
  container.className = 'matrix-input-container';
  container.style.gridTemplateColumns = `repeat(${size}, 1fr)`;
  
  for (let i = 0; i < size; i++) {
    for (let j = 0; j < size; j++) {
      const input = document.createElement('input');
      input.type = 'text';
      input.className = 'matrix-input-cell';
      input.placeholder = i === j ? '1' : '0';
      input.dataset.row = i;
      input.dataset.col = j;
      container.appendChild(input);
    }
  }
  
  matrixInput.appendChild(container);
}


// Circuit-based custom gate setup
function setupCircuitBuilder() {
  const numQubits = parseInt(circuitQubits.value);
  circuitGatesList.innerHTML = '<div class="small">No gates added to circuit yet.</div>';
  
  // Setup circuit gate selector
  setupCircuitGateSelector();
}

function setupCircuitGateSelector() {
  addCircuitGate.onclick = () => {
    const gateType = prompt('Enter gate type (X, Y, Z, H, S, T, CNOT, CZ, SWAP, CCNOT, etc.):');
    if (!gateType) return;

    const target = prompt(`Enter target qubit (0-${currentCustomGate.qubits-1}):`);
    if (target === null) return;

    const targetQubit = parseInt(target);
    if (isNaN(targetQubit) || targetQubit < 0 || targetQubit >= currentCustomGate.qubits) {
      alert('Invalid target qubit');
      return;
    }

    let gate = { type: gateType, params :[targetQubit] };

    // Handle controlled gates
    if (['CNOT', 'CZ', 'SWAP', 'CCNOT'].includes(gateType)) {
      let numControls = 1;
      if (gateType === 'CCNOT') numControls = 2; // Toffoli specifically
      if (gateType === 'SWAP') numControls = 0; // SWAP = two targets, special case

      for (let k = 0; k < numControls; k++) {
        const control = prompt(`Enter control qubit ${k+1} (0-${currentCustomGate.qubits-1}):`);
        if (control === null) return;

        const controlQubit = parseInt(control);
        if (isNaN(controlQubit) || controlQubit < 0 || controlQubit >= currentCustomGate.qubits) {
          alert('Invalid control qubit');
          return;
        }

        if (controlQubit === targetQubit) {
          alert('Control and target must be different');
          return;
        }
        if(controlQubit){
          if (gate.params.includes(controlQubit)) {
            alert('Duplicate control qubit');
            return;
          }
        }
        gate.params.push(controlQubit);
      }

      if (gateType === 'SWAP' ) {
        // SWAP needs two targets instead of control
        const swapTarget = prompt(`Enter second target qubit (0-${currentCustomGate.qubits-1}):`);
        if (swapTarget === null) return;

        const swapQubit = parseInt(swapTarget);
        if (isNaN(swapQubit) || swapQubit < 0 || swapQubit >= currentCustomGate.qubits || swapQubit === targetQubit) {
          alert('Invalid swap target qubit');
          return;
        }

        gate.params = swapQubit;
      }

    }

    // Handle rotation gates
    if (['Rx', 'Ry', 'Rz', 'Phase'].includes(gateType)) {
      const angle = prompt('Enter angle in degrees:');
      if (angle === null) return;

      const angleRad = parseFloat(angle) * Math.PI / 180;
      if (isNaN(angleRad)) {
        alert('Invalid angle');
        return;
      }

      gate.params = angleRad;
    }

    currentCustomGate.gates.push(gate);
    console.log(currentCustomGate);
    updateCircuitGatesList();
  };
}



function updateCircuitGatesList() {
  circuitGatesList.innerHTML = '';
  
  if (currentCustomGate.gates.length === 0) {
    circuitGatesList.innerHTML = '<div class="small">No gates added to circuit yet.</div>';
    return;
  }
  
  currentCustomGate.gates.forEach((gate, index) => {
    const item = document.createElement('div');
    item.className = 'gate-item';
    
    const left = document.createElement('div');
    left.className = 'gate-left';
    let desc = `${index + 1}. ${gate.type}`;
    if (gate.params?.length) {
      desc += ` (${gate.params.join(',')})`;
    }
    if (gate.angle !== undefined) {
      const deg = (gate.angle * 180 / Math.PI).toFixed(1);
      desc += `, ${deg}°`;
    }
    left.textContent = desc;
    
    const right = document.createElement('div');
    right.className = 'gate-right';
    
    const removeBtn = document.createElement('button');
    removeBtn.textContent = 'Remove';
    removeBtn.className = 'rm';
    removeBtn.onclick = () => {
      currentCustomGate.gates.splice(index, 1);
      updateCircuitGatesList();
    };
    
    right.appendChild(removeBtn);
    item.appendChild(left);
    item.appendChild(right);
    circuitGatesList.appendChild(item);
  });
}

// Control-based custom gate setup
function setupControlBuilder() {
  controlQubitsList.innerHTML = '';
  controlTarget.innerHTML = '';
  
  // Populate control qubits (checkboxes)
  for (let i = 0; i < nQ; i++) {
    const item = document.createElement('div');
    item.className = 'control-qubit-item';
    
    const checkbox = document.createElement('input');
    checkbox.type = 'checkbox';
    checkbox.id = `control_${i}`;
    checkbox.value = i;
    
    const label = document.createElement('label');
    label.htmlFor = `control_${i}`;
    label.textContent = `q${i}`;
    
    item.appendChild(checkbox);
    item.appendChild(label);
    controlQubitsList.appendChild(item);
  }
  
  // Populate target qubit
  for (let i = 0; i < nQ; i++) {
    const option = document.createElement('option');
    option.value = i;
    option.textContent = `q${i}`;
    controlTarget.appendChild(option);
  }
  
  // Show/hide angle input for rotation gates
  if (['Rx', 'Ry', 'Rz', 'Phase'].includes(baseGate.value)) {
    controlAngleDiv.classList.remove('hidden');
  } else {
    controlAngleDiv.classList.add('hidden');
  }
}

// Matrix validation
function validateCustomMatrix() {
  const name = customGateName.value.trim();
  if (!name) {
    showValidation(matrixValidation, 'Please enter a gate name', 'error');
    return false;
  }
  
  const size = parseInt(matrixSize.value);
  if (isNaN(size) || size < 1 || size > 8) {
    showValidation(matrixValidation, 'Invalid matrix size', 'error');
    return false;
  }
  if(size===1){

  }
  const matrix = [];
  const inputs = matrixInput.querySelectorAll('.matrix-input-cell');
  
  if (inputs.length !== size * size) {
    showValidation(matrixValidation, 'Matrix input not properly initialized', 'error');
    return false;
  }
  
  for (let i = 0; i < size; i++) {
    matrix[i] = [];
    for (let j = 0; j < size; j++) {
      const input = matrixInput.querySelector(`[data-row="${i}"][data-col="${j}"]`);
      if (!input) {
        showValidation(matrixValidation, `Missing input at position (${i},${j})`, 'error');
        return false;
      }
      
      const value = input.value.trim() || (i === j ? '1' : '0');
      
      try {
        matrix[i][j] = parseComplex(value);
      } catch (e) {
        showValidation(matrixValidation, `Invalid complex number at position (${i},${j}): "${value}". Use format like: 1, 0, i, 1+i, 0.5-0.3i`, 'error');
        return false;
      }
    }
  }
  
  // Check if matrix is unitary
  try {
    if (!isUnitary(matrix)) {
      showValidation(matrixValidation, 'Matrix is not unitary (U†U ≠ I). Please ensure the matrix represents a valid quantum gate.', 'error');
      return false;
    }
  } catch (e) {
    showValidation(matrixValidation, `Error validating matrix: ${e.message}`, 'error');
    return false;
  }
  const target = parseInt(controlTarget.value);
  // Save custom gate
  customGates[name] = {
    type: 'CUSTOM',
    name: name,
    customType : 'CUSTOM_MATRIX',
    matrix: matrix,
    size: size,
    params : [target]
  };
  console.log(customGates[name]);
  showValidation(matrixValidation, `Custom gate "${name}" created successfully!`, 'success');
  return true;
}

//Circuit validation
function validateCustomCircuit() {
  const name = circuitGateName.value.trim();
  if (!name) {
    showValidation(circuitValidation, 'Please enter a gate name', 'error');
    return false;
  }
  
  if (currentCustomGate.gates.length === 0) {
    showValidation(circuitValidation, 'Please add at least one gate to the circuit', 'error');
    return false;
  }
  const nqubits = circuitQubits.value.trim();
  let num = prompt("enter the qubits numbers on which these circuit should applied ");
  let params = num.split(" ").map(Number);
  // Save custom gate
  customGates[name] = {
    type: 'CUSTOM',
    name: name,
    customType : 'CUSTOM_CIRCUIT',
    qubits: currentCustomGate.qubits,
    subGates: [...currentCustomGate.gates],
    params : params
  };
  console.log(customGates[name]);
  
  showValidation(circuitValidation, `Custom circuit "${name}" created successfully!`, 'success');
  return true;
}

// Control gate validation
function validateCustomControl() {
  const name = controlGateName.value.trim();
  if (!name) {
    showValidation(controlValidation, 'Please enter a gate name', 'error');
    return false;
  }
  
  const selectedControls = Array.from(controlQubitsList.querySelectorAll('input[type="checkbox"]:checked'))
    .map(cb => parseInt(cb.value));
  
  if (selectedControls.length === 0) {
    showValidation(controlValidation, 'Please select at least one control qubit', 'error');
    return false;
  }
  
  const target = parseInt(controlTarget1.value);
  if (selectedControls.includes(target)) {
    showValidation(controlValidation, 'Target qubit cannot be the same as control qubits', 'error');
    return false;
  }
  
  const type = baseGate.value;

  let gate = { type , params: []};

  const angle = ['Rx', 'Ry', 'Rz', 'Phase'].includes(type) ? parseFloat(controlAngle.value || 0) * Math.PI / 180 : null;
  if (['X','Y','Z','H','S','Sdg','T','Tdg','Rx','Ry','Rz','Phase','MEASURE'].includes(type)){
    gate.params = [target];
    if (['Rx','Ry','Rz','Phase'].includes(type)){
      gate.angle = angle; // store radians
    }
  }
  // Save custom gate
  customGates[name] = {
    type: 'CUSTOM',
    name: name,
    customType : 'CUSTOM_CONTROL',
    subGates: [gate],
    params : [...selectedControls, target],
    angle: angle
  };

  console.log(customGates[name]);  
  showValidation(controlValidation, `Custom control gate "${name}" created successfully!`, 'success');
  return true;
}

// Helper functions
function parseComplex(str) {
  if (!str || typeof str !== 'string') return c(0, 0);
  
  str = str.trim();
  if (str === '0' || str === '') return c(0, 0);
  if (str === '1') return c(1, 0);
  if (str === 'i') return c(0, 1);
  if (str === '-i') return c(0, -1);
  
  // Handle simple real numbers
  if (!str.includes('i') && !str.includes('j')) {
    const real = parseFloat(str);
    if (!isNaN(real)) return c(real, 0);
  }
  
  // Parse a+bi format (more flexible)
  const match = str.match(/^([+-]?\d*\.?\d*)\s*([+-]?\d*\.?\d*)?[ij]?$/);
  if (!match) {
    // Try alternative parsing
    const altMatch = str.match(/^([+-]?\d*\.?\d*)\s*([+-]?\d*\.?\d*)?[ij]$/);
    if (altMatch) {
      const real = parseFloat(altMatch[1]) || 0;
      const imag = parseFloat(altMatch[2]) || 1;
      return c(real, imag);
    }
    throw new Error(`Invalid complex number format: "${str}"`);
  }
  
  const real = parseFloat(match[1]) || 0;
  const imag = parseFloat(match[2]) || 0;
  
  return c(real, imag);
}

function isUnitary(matrix) {
  const n = matrix.length;
  const identity = Array(n).fill().map(() => Array(n).fill(c(0, 0)));
  for (let i = 0; i < n; i++) identity[i][i] = c(1, 0);
  
  // Compute U†U
  const product = Array(n).fill().map(() => Array(n).fill(c(0, 0)));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      for (let k = 0; k < n; k++) {
        // Ensure matrix values are valid complex numbers
        const conjVal = math.conj(matrix[k][i] || c(0, 0));
        const matrixVal = matrix[k][j] || c(0, 0);
        const multResult = math.multiply(conjVal, matrixVal);
        product[i][j] = math.add(product[i][j] || c(0, 0), multResult);
      }
    }
  }
  
  // Check if product equals identity (within tolerance)
  const tolerance = 1e-10;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const diff = math.subtract(product[i][j] || c(0, 0), identity[i][j] || c(0, 0));
      if (math.abs(diff) > tolerance) return false;
    }
  }
  return true;
}

function showValidation(element, message, type) {
  element.innerHTML = `<div class="validation-message validation-${type}">${message}</div>`;
}


// Event listeners for custom gates
matrixSize.addEventListener('change', setupMatrixInput);
validateMatrix.addEventListener('click', validateCustomMatrix);
circuitQubits.addEventListener('change', setupCircuitBuilder);
validateCircuit.addEventListener('click', validateCustomCircuit);
baseGate.addEventListener('change', () => {
  if (['Rx', 'Ry', 'Rz', 'Phase'].includes(baseGate.value)) {
    controlAngleDiv.classList.remove('hidden');
  } else {
    controlAngleDiv.classList.add('hidden');
  }
});
validateControl.addEventListener('click', validateCustomControl);

// ---------- Interactivity Enhancements ----------
function toggleTheme(){
  const isDark = document.body.classList.toggle('dark');
  localStorage.setItem('theme', isDark ? 'dark' : 'light');
  if (themeToggleBtn) themeToggleBtn.textContent = isDark ? '☀️' : '🌙';
  // update plot backgrounds immediately
  refreshPlotBackgrounds();
}
function initThemeFromStorage(){
  const saved = localStorage.getItem('theme');
  const prefersDark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
  const useDark = saved ? (saved === 'dark') : prefersDark;
  if (useDark) document.body.classList.add('dark');
  if (themeToggleBtn) themeToggleBtn.textContent = useDark ? '☀️' : '🌙';
}

function getPlotColors(){
  const isDark = document.body.classList.contains('dark');
  return isDark
    ? { paper:'#8ca3d3', plot:'#8ca3d3', scene:'#8ca3d3' }
    : { paper:'rgb(246, 246, 248)', plot:'rgb(247, 247, 247)', scene:'rgb(249, 249, 250)' };
}

function refreshPlotBackgrounds(){
  const colors = getPlotColors();
  // Q-sphere
  const qDiv = document.getElementById('qsphereDiv');
  if (qDiv && qDiv.data) {
    try { Plotly.relayout(qDiv, { paper_bgcolor: colors.paper, plot_bgcolor: colors.plot, 'scene.bgcolor': colors.scene }); } catch(e) {}
  }
  // Bloch spheres
  document.querySelectorAll('.bloch-canvas').forEach((el)=>{
    try { Plotly.relayout(el, { paper_bgcolor: colors.paper, plot_bgcolor: colors.plot, 'scene.bgcolor': colors.scene }); } catch(e) {}
  });
}

function setControlsDisabled(disabled){
  const controls = document.querySelectorAll('button, input, select');
  controls.forEach(el => {
    if (el.id === 'themeToggle') return; // keep theme usable
    el.disabled = disabled;
    el.style.opacity = disabled ? 0.7 : 1;
    el.style.pointerEvents = disabled ? 'none' : 'auto';
  });
}
function showLoading(show){
  if (!loadingSpinner) return;
  loadingSpinner.style.display = show ? 'block' : 'none';
}
async function onBackendRunClickWrap(ev){
  // mimic original click handler but add loading/disable UX
  try{
    if (isRunningBackend) return; // prevent double runs
    isRunningBackend = true;
    showLoading(true);
    setControlsDisabled(true);
    // trigger original listener logic by dispatching a new click is redundant; call inline logic instead
    await backendRunCore();
  } finally {
    showLoading(false);
    setControlsDisabled(false);
    isRunningBackend = false;
  }
}
async function backendRunCore(){
  // Factor out body of existing handler to reuse with UX wrapper
  const nQ_local = nQ;
  if(nQ_local<6){
    const payload = { numQubits: nQ_local, gates: gateSequence, initialStates: initStates };
    const res = await fetch("https://qsv-3xax.onrender.com/run", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload)
    });
    console.log("👉 Sending to backend:", payload);
    const data = await res.json();
    document.getElementById("backendResults").textContent =
      "Backend Counts:\n" + JSON.stringify(data.counts, null, 2) +
      "\n\nQASM:\n" + data.qasm;
    drawHistogram(data.counts);
  } else {
    alert("we used single qubit Quantum State tomography for qubits more than 5. The results may have some eroor due to little noise ");
    targetQubit = parseInt(document.getElementById("qubitSelect").value);
    if (isNaN(targetQubit)) {
      alert("Please select a valid target qubit!");
      return;
    }
    const payload = { numQubits: nQ_local, gates: gateSequence, initialStates: initStates, targetQubit };
    const res = await fetch("https://qsv-3xax.onrender.com/run", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload)
    });
    const data = await res.json();
    document.getElementById("backendResults").textContent =
      `Bloch : (${data.blochs.x}, ${data.blochs.y}, ${data.blochs.z})`+
      "\n\nRhos:\n" + formatComplexMatrix(data.rho);
    MathJax.typeset();
    const wrapper = document.createElement('div');
    wrapper.className = 'bloch-wrapper';
    const block = document.createElement('div');
    block.id = 'bloch_' + targetQubit;
    block.className = 'bloch-canvas';
    wrapper.appendChild(block);
    const props = document.createElement('div');
    props.className = 'bloch-properties';
    const x = Number(data.blochs.x.toFixed(3));
    const y = Number(data.blochs.y.toFixed(3));
    const z = Number(data.blochs.z.toFixed(3));
    const e = cleanFixed(qubitEntropy(x,y,z).toFixed(3));
    props.innerHTML = `
      <h3>Qubit ${targetQubit}</h3>
      <p>Bloch vector:(${x}, ${y}, ${z})</p>
      <p>Entropy(mixedness): ${e}</p>
      <p>Purity: ${(1 + x * x + y * y + z * z) / 2}</p>
      <p>Measurement probabilities(|0>,|1>): 
         ${formatReal(data.rho[0][0], 3)}, 
         ${formatReal(data.rho[1][1], 3)}
      </p>
    `;
    wrapper.appendChild(props);
    blochSpheresDiv.appendChild(wrapper);
    const r = Math.sqrt(x**2 + y**2 + z**2);
    if(r< 1e-6){
      plotBloch(block.id, {x:0,y:0,z:0}, targetQubit,true);
    }
    plotBloch(block.id,{x,y,z},targetQubit,true);
  }
}

// Keyboard shortcuts
window.addEventListener('keydown', (e)=>{
  if (e.target && (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT' || e.target.isContentEditable)) return;
  if (e.key === 'a' || e.key === 'A') { document.getElementById('btnAddGate').click(); }
  if (e.key === 'u' || e.key === 'U') { btnUndo.click(); }
  if (e.key === 'c' || e.key === 'C') { btnClearGates.click(); }
  if (e.key === 'r' || e.key === 'R') {
    if (e.shiftKey) { runCircuitBtn.click(); }
    else { btnRun.click(); }
  }
});
   



