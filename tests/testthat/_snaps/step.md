# num_comp is reduced below the number of columns

    Code
      trained <- prep(step_manifold_sculpting(rec, x, y, neighbors = 14, num_comp = 3,
        seed = 1, options = list(max_iter = 20)))
    Condition
      Warning:
      `num_comp` was reduced to 1, one less than the number of selected columns.

# missing columns at bake time are an error

    Code
      bake(trained, new_data = test[c("x", "y", "label")])
    Condition
      Error in `bake()`:
      x The following required columns are missing from `new_data`: `z`.
      i These columns have one of the following roles, which are required at `bake()` time: `predictor`.

# invalid arguments are errors

    Code
      prep(step_manifold_sculpting(rec, x, y, z, options = list(foo = 1)))
    Condition
      Error in `step_manifold_sculpting()`:
      Caused by error in `prep()`:
      ! `options` has unknown element "foo".
      i Allowed elements are "nstart", "max_iter", "patience", "nn_method", and "verbose".
    Code
      prep(step_manifold_sculpting(rec, x, y, z, squishing_rate = 2))
    Condition
      Error in `step_manifold_sculpting()`:
      Caused by error in `prep()`:
      ! `squishing_rate` must be a single number in (0, 1), not 2.
    Code
      prep(step_manifold_sculpting(rec, x, y, z, neighbors = -1))
    Condition
      Error in `step_manifold_sculpting()`:
      Caused by error in `prep()`:
      ! `neighbors` must be a single whole number >= 1, not -1.
    Code
      prep(step_manifold_sculpting(rec, x, y, z, prefix = 1))
    Condition
      Error in `step_manifold_sculpting()`:
      Caused by error in `prep()`:
      ! `prefix` must be a single string, not 1.
    Code
      prep(step_manifold_sculpting(rec, label))
    Condition
      Error in `step_manifold_sculpting()`:
      Caused by error in `prep()`:
      x All columns selected for the step should be double or integer.
      * 1 factor variable found: `label`

# printing

    Code
      print(step)
    Message
      
      -- Recipe ----------------------------------------------------------------------
      
      -- Inputs 
      Number of variables by role
      outcome:   1
      predictor: 3
      
      -- Operations 
      * Manifold Sculpting extraction with: x, y, z

---

    Code
      prep(step)
    Message
      
      -- Recipe ----------------------------------------------------------------------
      
      -- Inputs 
      Number of variables by role
      outcome:   1
      predictor: 3
      
      -- Training information 
      Training data contained 300 data points and no incomplete rows.
      
      -- Operations 
      * Manifold Sculpting extraction with: x, y, z | Trained

---

    Code
      step_manifold_sculpting(rec, x, y, z, num_comp = 0)
    Message
      
      -- Recipe ----------------------------------------------------------------------
      
      -- Inputs 
      Number of variables by role
      outcome:   1
      predictor: 3
      
      -- Operations 
      * Manifold Sculpting was not conducted for: x, y, z

