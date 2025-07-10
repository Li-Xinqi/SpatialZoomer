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
  '<strong>Dataset 1:</strong> Xenium In Situ, Alzheimer\'s Disease Mouse Model<br><strong>Description:</strong> This dataset comprises six FFPE tissue sections, organized into two columns representing wild-type (WT) and CRND8 APP-overexpressing (TgCRND8) transgenic male mouse brains. The three rows correspond to time points of 2.5, 5.7, and 13.2 months of age for WT mice, and 2.5, 5.7, and 17.9 months for TgCRND8 mice. In the transgenic group, these time points reflect progressive stages of Alzheimer’s disease pathology, specifically mild, moderate, and advanced amyloid-β (Aβ) deposition.<br><strong>Cells:</strong> 349,063 &nbsp;&nbsp; <strong>Genes:</strong> Xenium Mouse Brain Gene Expression Panel plus 99 Custom Genes (347 genes)<br><strong>Original Data Download:</strong> <a href="https://www.10xgenomics.com/datasets/xenium-in-situ-analysis-of-alzheimers-disease-mouse-model-brain-coronal-sections-from-one-hemisphere-over-a-time-course-1-standard" target="_blank">Link</a>',
  "Vizgen MERFISH, :.",         // Dataset2
  "Xenium Prime.",       // Dataset3
  "Xenium Prime, ."              // Dataset4
];

// 二维参数数组，按 dataset 索引存储对应的 scale 和 resolution 列表
const scaleOptionsList = [
  [0.0, 1.1, 3.0, 4.5, 6.0, 8.0, 10.5, 14.0, 25.0, 40.0],       // Dataset1 scales
  [0.0, 1, 3, 4, 5.5, 7.5, 9.5, 13.5, 25.0, 40.0],         // Dataset2 scales
  [0.0, 0.6, 1.2, 1.6, 3.0, 4.0, 5.0, 6.5, 8.5, 13.0, 25.0, 40.0],       // Dataset3 scales
  [0.0, 1.0, 2.5, 3.5, 4.5, 6.0, 8.0, 13.0, 25.0, 40.0]              // Dataset4 scales
];

const resolutionOptionsList = [
  ['1', '1.5', '2'],          // Dataset1 resolutions
  ['1', '0.8', '1.2'],               // Dataset2 resolutions
  ['0.8', '1', '0.1'],     // Dataset3 resolutions
  ['0.8', '1', '0.6']           // Dataset4 resolutions
];

const resolutionOptionsValue = [
  ['1', '1.5', '2'],          // Dataset1 resolutions
  ['1', '0.8', '1.2'],               // Dataset2 resolutions
  ['0.8', '1', '0.1'],     // Dataset3 resolutions
  ['0.8', '1', '0.6']           // Dataset4 resolutions
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
  scaleValue.textContent = scaleOptionsList[currentDataset-1][idx].toFixed(1);
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
  scaleValue.textContent = scales[0].toFixed(1);
}

function updateDemo() {
  if (!currentDataset) return;
  desc.innerHTML = datasetDescriptionList[currentDataset-1];
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
  img3.src = `${baseurl}/assets/images/data/demo${currentDataset}_scale_${scale}_reso${res}.jpg`;
}