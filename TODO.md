#### Parallel Tempering

- `thermalize()` should be able to use `Freezer` and `ConstantT` (implement abstract type)
- `parallel_tempering!()` should perform `init_edges!()` when necessary
- should use `ConstantT` by default?
- dynamic temperature range updates
  - https://arxiv.org/pdf/1501.05823.pdf
  - CMBP script


#### General

- move to using my new Lattice library (it's nicer and more straight forward)
- create files after the simulation... 


#### data collection

- this should really be split from `measure!()`, it creates a lot of unnecessary clutter
