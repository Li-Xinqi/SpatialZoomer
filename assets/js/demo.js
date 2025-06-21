---
---
// const baseurl = window.jekyllBaseurl || '';
var baseurl = "{{ site.baseurl }}";
const buttons = document.querySelectorAll('.dataset-buttons button');
const desc = document.getElementById('dataset-description');
const img1 = document.getElementById('img1');
const img2 = document.getElementById('img2');
const img3 = document.getElementById('img3');
const scaleRange = document.getElementById('scaleRange');
const scaleValue = document.getElementById('scaleValue');
const scaleOptions = [0.1, 0.2, 0.5, 1, 2, 5, 10];
const resolutionSelect = document.getElementById('resolutionSelect');
const resolutionOptions = ['1', '2', '3'];
let currentDataset = null;

document.addEventListener('DOMContentLoaded', () => {
  currentDataset = buttons[0].dataset.dataset; // 默认选择第一个数据集
  updateDemo();
});

buttons.forEach(btn => btn.addEventListener('click', () => {
  currentDataset = btn.dataset.dataset;
  updateDemo();
}));

scaleRange.addEventListener('input', () => {
  const idx = parseInt(scaleRange.value, 10);
  scaleValue.textContent = scaleOptions[idx];
  updateImages();
});

// 页面加载时动态填充 resolution 选项，并清空旧选项以防重复
function populateResolutionOptions() {
  resolutionSelect.innerHTML = '';
  resolutionOptions.forEach(res => {
    const opt = document.createElement('option');
    opt.value = res;
    opt.textContent = res;
    resolutionSelect.appendChild(opt);
  });
}

resolutionSelect.addEventListener('change', updateImages);

function updateDemo() {
  if (!currentDataset) return;
  desc.textContent = `Dataset ${currentDataset} description.`;
  img1.src = `${baseurl}/assets/images/demo${currentDataset}_1.png`;
  img2.src = `${baseurl}/assets/images/demo${currentDataset}_2.png`;
  scaleRange.value = 0;
  scaleValue.textContent = scaleOptions[scaleRange.value];
  // 动态填充 resolution 下拉选项
  populateResolutionOptions()
  resolutionSelect.value = 1;
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