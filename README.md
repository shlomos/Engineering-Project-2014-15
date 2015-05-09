<h1>3D Anatomy Segmentation using Perceptual Computing</h1>

<p>The delineation and segmentation of anatomical structures and pathologies in volumetric scans is a challenging and time-consuming task that significantly hampers the use of 3D models in the clinic.
We present a semi-automatic tool for 3D segmentation in volumetric medical scans using natural input from the user. The input consists of hand motions and gestures acquired using the Leap Motion Device. The input is a 3D approximate scribble across several slices of the CT scan; the output consists of both segmentation of the corrected model and the 3D mesh of it. </p>
<p>Our method consists of four steps: 1) acquire two initial approximate 2D scribbles of the anatomy to be segmented and of the background; 2) online segmentation of the structure of interest with the Grow-Cut segmentation algorithm over the given initial scribble and real-time 2D correction; 3) 3D correction of the segmentation leaks over the mesh using natural input, and; 4) iteratively repeat stages 2-3 until a satisfactory result is obtained. </p>
<p>To evaluate our method, we performed both a qualitative analysis on healthy soft-tissue anatomies and quantitative examination on scans of liver tumors from 10 patients. The experimental results for the quantitative examination for liver tumors yield <b>Volume Overlap Error of 15.1% (std : 4.59%)</b> and  <b>Volume Similarity Error of 8.74% (std: 4.77%)</b>.</p>

<img src='Group209 3D Anatomy Segmentation using Perceptual Computing.pic1.JPG'>

<p><strong>
<h4>Project members</h4>
<br/>
Shlomo    Shenzis		shlomo.shenzis@mail.huji.ac.il<br/>
Moshe    Samson		moshe.samson@mail.huji.ac.il
<br/><br/>
<h4>Supervisors</h4>
<br/>
Prof. Leo Joskowicz, The Hebrew University of Jerusalem<br/>
Mr. Rafael Vivanti, The Hebrew University of Jerusalem 
</strong></p>
