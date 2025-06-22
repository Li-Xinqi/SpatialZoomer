const baseurl = window.jekyll_baseurl || "";
const buttons = document.querySelectorAll('.dataset-buttons button');
const desc = document.getElementById('dataset-description');
const img1 = document.getElementById('img1');
const img2 = document.getElementById('img2');
const img3 = document.getElementById('img3');
const scaleRange = document.getElementById('scaleRange');
const scaleValue = document.getElementById('scaleValue');
const resolutionSelect = document.getElementById('resolutionSelect');

// 二维参数数组，按 dataset 索引存储对应的 scale 和 resolution 列表
const datasetDescriptionList = [
  "Dataset 1 description.",       // Dataset1
  "Dataset 2 description.",         // Dataset2
  "Dataset 3 description.",       // Dataset3
  "Dataset 4 description."              // Dataset4
];

// 二维参数数组，按 dataset 索引存储对应的 scale 和 resolution 列表
const scaleOptionsList = [
  [0.1, 0.2, 0.5, 1, 2, 5, 10],       // Dataset1 scales
  [0.2, 0.5, 1, 2],         // Dataset2 scales
  [0.5, 1, 2, 5, 10],       // Dataset3 scales
  [1, 2, 4, 8]              // Dataset4 scales
];

const resolutionOptionsList = [
  ['1', '2', '3'],          // Dataset1 resolutions
  ['a', 'b'],               // Dataset2 resolutions
  ['x', 'y', 'z', 'w'],     // Dataset3 resolutions
  ['low', 'high']           // Dataset4 resolutions
];

const resolutionOptionsValue = [
  ['1', '2', '3'],    // Dataset1 actual values
  ['10', '20'],             // Dataset2
  ['100', '200', '300', '400'], // Dataset3
  ['1.0', '2.0']            // Dataset4
];

let currentDataset = null;

document.addEventListener('DOMContentLoaded', () => {
  currentDataset = parseInt(buttons[0].dataset.dataset); // 默认选择第一个数据集
  updateDemo();
});

buttons.forEach(btn => btn.addEventListener('click', () => {
  currentDataset = parseInt(btn.dataset.dataset);
  updateDemo();
}));

scaleRange.addEventListener('input', () => {
  const idx = parseInt(scaleRange.value, 10);
  scaleValue.textContent = scaleOptionsList[currentDataset-1][idx];
  updateImages();
});

resolutionSelect.addEventListener('change', updateImages);

// 填充 resolution 下拉选项，清空旧选项以防重复
function populateResolutionOptions(texts, values) {
  resolutionSelect.innerHTML = '';
  texts.forEach((text, idx) => {
    const opt = document.createElement('option');
    opt.value = values[idx];
    opt.textContent = text;
    resolutionSelect.appendChild(opt);
  });
  resolutionSelect.value = texts[0]; // 默认选择第一个选项
}

// 初始化 scaleRange
function initializeScales(scales) {
  scaleRange.min = 0;
  scaleRange.max = scales.length - 1;
  scaleRange.step = 1;
  scaleRange.value = 0;
  scaleValue.textContent = scales[0];
}

function updateDemo() {
  if (!currentDataset) return;
  desc.textContent = datasetDescriptionList[currentDataset-1];
  img1.src = `${baseurl}/assets/images/demo${currentDataset}_1.png`;
  img2.src = `${baseurl}/assets/images/demo${currentDataset}_2.png`;
  // 初始化 scaleRange
  initializeScales(scaleOptionsList[currentDataset-1]);
  // 动态填充 resolution 下拉选项
  populateResolutionOptions(resolutionOptionsList[currentDataset-1], resolutionOptionsValue[currentDataset-1]);
  updateImages();
}

// function updateScale() {
//   if (!currentDataset) return;
//   const scale = scaleRange.value;
//   img3.src = `/assets/images/demo${currentDataset}_3_scale${scale}_reso${reso}.png`;
// }

function updateImages() {
  if (!currentDataset) return;
  const scale = scaleValue.textContent;
  const res = resolutionSelect.value;
  img3.src = `${baseurl}/assets/images/data/demo${currentDataset}_scale${scale}_reso${res}.png`;
}