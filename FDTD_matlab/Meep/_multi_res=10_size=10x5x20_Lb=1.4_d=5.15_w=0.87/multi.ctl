(set! resolution 10)
(set! force-complex-fields? true)
(set! pml-layers (list
   (make pml (thickness 1))
))
(set! sources (list
   (make source (src (make continuous-src (frequency 0.26316)))(component Ex)(center 0 0 8.05)(size 10 5 0))
   (make source (src (make continuous-src (frequency 0.26316)))(component Ey)(center 0 0 9)(size 10 5 0))
))
(set! geometry-lattice (make lattice (size 10 5 20)))
(set! geometry (list
   (make block (material (make dielectric (index 2.3812)))(center -4 -1.8 0.575)(size 2 0.83 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center -4 -0.4 0.575)(size 2 0.83 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center -4 1 0.575)(size 2 0.83 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center -4 2.4 0.575)(size 2 0.83 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center -2 -1.8 0.575)(size 2 0.85 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center -2 -0.4 0.575)(size 2 0.85 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center -2 1 0.575)(size 2 0.85 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center -2 2.4 0.575)(size 2 0.85 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 0 -1.8 0.575)(size 2 0.87 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 0 -0.4 0.575)(size 2 0.87 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 0 1 0.575)(size 2 0.87 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 0 2.4 0.575)(size 2 0.87 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 2 -1.8 0.575)(size 2 0.89 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 2 -0.4 0.575)(size 2 0.89 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 2 1 0.575)(size 2 0.89 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 2 2.4 0.575)(size 2 0.89 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 4 -1.8 0.575)(size 2 0.91 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 4 -0.4 0.575)(size 2 0.91 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 4 1 0.575)(size 2 0.91 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 4 2.4 0.575)(size 2 0.91 5.15)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
   (make block (material (make dielectric (index 2.3812)))(center 0.575 0 -16)(size 20 10 28)(e1 1 0 0)(e2 0 1 0)(e3 0 0 1))
))

(run-until 200
   (at-beginning
      output-epsilon
   )
   (at-end
      (in-volume 
         (volume (center 0 0 -2) (size 10 5 0))
         (to-appended "Ex_z=0" output-efield-x)
         (to-appended "Ey_z=0" output-efield-y)
      )
   )
   (at-end
      (in-volume 
         (volume (center 0 0 -10) (size 10 5 0))
         (to-appended "Ex_z=8" output-efield-x)
         (to-appended "Ey_z=8" output-efield-y)
      )
   )
)
