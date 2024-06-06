# Stabler Neo-Hookean Simulation: <br> Absolute Eigenvalue Filtering for Projected Newton

### [Project Page](https://www.cs.columbia.edu/cg/local-deformation/)  | [Paper](http://www.cs.columbia.edu/cg/abs-eval/paper_low_res.pdf)

<img src="https://github.com/honglin-c/abs-eval/blob/main/.github/images/teaser.png" width="800">

Official implementation for the paper:
> **[Stabler Neo-Hookean Simulation: Absolute Eigenvalue Filtering for Projected Newton](https://www.cs.columbia.edu/cg/local-deformation/)**  
> [Honglin Chen](https://www.cs.columbia.edu/~honglinchen/)<sup>1</sup>, 
[Hsueh-Ti Derek Liu](https://www.dgp.toronto.edu/~hsuehtil/)<sup>3,</sup><sup>4</sup>, 
[David I.W. Levin](http://www.cs.toronto.edu/~diwlevin/)<sup>2,</sup><sup>5</sup>, 
[Changxi Zheng](http://www.cs.columbia.edu/~cxz/)<sup>1</sup>, 
[Alec Jacobson](https://www.cs.toronto.edu/~jacobson/)<sup>2,</sup><sup>6</sup>
<br>
> <sup>1</sup>Columbia University, 
<sup>2</sup>University of Toronto,  
<sup>3</sup>Roblox Research, 
<sup>4</sup>University of British Columbia, 
<br> 
&nbsp; &nbsp; <sup>5</sup>NVIDIA,
<sup>6</sup>Adobe Research
<br>
> SIGGRAPH 2024 (Conference Track)


## Installation
To build the code, please run
```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
```

## Experiments

To run the code with our absolute eigenvalue projection strategy, please run, e.g.,
```
sh ../scripts/micky_teaser/abs.sh
```
We provide several examples in `scripts/`. 

### Comparisons
To run the same example with the eigenvalue clamping strategy, run `clamp.sh` under the same folder (or add `--clamp` option after the command in `abs.sh`):
```
sh ../scripts/micky_teaser/clamp.sh
```
The default eigenvalue clamping algorithm uses 0 as the clamping threshold. To use a different clamping threshold (e.g., a small positive number), add `--epsilon [threshold]` option after the command. 

## Optional arguments

For more options, please see
```
./example --help
```

<details>
<summary>
<h2>What do we modify in TinyAD to achieve absolute eigenvalue projection? </h2>
</summary>

We clone TinyAD (https://github.com/patr-schm/TinyAD/blob/29417031c185b6dc27b6d4b684550d844459b735D) to the project folder,
then comment out and change [lines 71-75 in `TinyAD/include/TinyAD/Utils/HessianProjection.hh`]([https://github.com/patr-schm/TinyAD/blob/29417031c185b6dc27b6d4b684550d844459b735/include/TinyAD/Utils/HessianProjection.hh#L71-L75]) to:
```
  if (_eigenvalue_eps < 0) {
      // project to absolute value if the eigenvalue threshold is less than 0
      if (D(i, i) < 0)
      {
          D(i, i) = -D(i, i);
          all_positive = false;
      }
  }
  else {
      // project to epsilon otherwise
      if (D(i, i) < _eigenvalue_eps)
      {
          D(i, i) = _eigenvalue_eps;
          all_positive = false;
      }
  }
```
So we project to absolute value if the eigenvalue threshold is less than 0, and to a small value epsilon (e.g., 0 or 1e-3) otherwise.

</details>

<details>
<summary>
<h2>Yet another way to verify our method (which is slightly more complicated) </h2>
</summary>

1. Git clone [Hobak](https://github.com/theodorekim/HOBAKv1/blob/8420c51b795735d8fb912e0f8810f935d96fb636) in a different folder.

2. Change [lines 139-141](https://github.com/theodorekim/HOBAKv1/blob/8420c51b795735d8fb912e0f8810f935d96fb636/src/Hyperelastic/Volume/SNH.cpp#L139-L141) in `src/Hyperelastic/Volume/SNH.cpp` to:
```
  for (int i = 0; i < 9; i++)
      if (eigenvalues(i) < 0.0)
          eigenvalues(i) = -eigenvalues(i);  
```

3. Change line 54 in `src/Scenes/QUASISTATIC_STRETCH.h` to:
```
  _hyperelastic = new VOLUME::SNH(1.0, 100.0);
```
So we increase the Poisson's ratio to roughly 0.495.

3. Change lines 80-84 in `src/Scenes/QUASISTATIC_STRETCH.h` to:
```
  if (_frameNumber < 1)
  {
    _kinematicShapes[0]->translation()[2] -= 0.5;
    return;
  }
```
So we create a large initial deformation (feel free to try an even larger deformation).

4. Uncomment [line 287](https://github.com/theodorekim/HOBAKv1/blob/8420c51b795735d8fb912e0f8810f935d96fb636/projects/simulateScene/simulateScene.cpp#L287) in `projects/simulateScene/simulateScene.cpp` to choose the quasistatic test:
```
  scene = new QUASISTATIC_STRETCH();
```

5. Run `make mac` or `make linux` to run the test (see HOBAK's readme).

* Note: HOBAK doesn't include a line search but our method seems to often work fine without one.

</details>

## Citation
```
@inproceedings{chen2024abs_eval,
      title={Stabler Neo-Hookean Simulation: Absolute Eigenvalue Filtering for Projected Newton},
      author={Honglin Chen and Hsueh-Ti Derek Liu and David I.W. Levin and Changxi Zheng and Alec Jacobson},
      booktitle = {ACM SIGGRAPH 2024 Conference Proceedings},
      year = {2024}
  }
```