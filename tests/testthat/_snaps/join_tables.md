# Joining tables from nt_descriptive by group and overall

    Code
      as.data.frame(tab)
    Output
                        Variable Species: setosa (n = 50) Species: virginica (n = 50)
      1             Sepal.Length                                                     
      2             \t Mean ± SD              5.01 ± 0.35                 6.59 ± 0.64
      3  \t Median (Q25% ; Q75%)       5.00 (4.80 ; 5.20)          6.50 (6.23 ; 6.90)
      4    \t Median (Min ; Max)       5.00 (4.30 ; 5.80)          6.50 (4.90 ; 7.90)
      5               \t Missing                        0                           0
      6              Sepal.Width                                                     
      7             \t Mean ± SD              3.43 ± 0.38                 2.97 ± 0.32
      8  \t Median (Q25% ; Q75%)       3.40 (3.20 ; 3.68)          3.00 (2.80 ; 3.18)
      9    \t Median (Min ; Max)       3.40 (2.30 ; 4.40)          3.00 (2.20 ; 3.80)
      10              \t Missing                        0                           0
      11            Petal.Length                                                     
      12            \t Mean ± SD              1.46 ± 0.17                 5.55 ± 0.55
      13 \t Median (Q25% ; Q75%)       1.50 (1.40 ; 1.58)          5.55 (5.10 ; 5.88)
      14   \t Median (Min ; Max)       1.50 (1.00 ; 1.90)          5.55 (4.50 ; 6.90)
      15              \t Missing                        0                           0
      16             Petal.Width                                                     
      17            \t Mean ± SD              0.25 ± 0.11                 2.03 ± 0.27
      18 \t Median (Q25% ; Q75%)       0.20 (0.20 ; 0.30)          2.00 (1.80 ; 2.30)
      19   \t Median (Min ; Max)       0.20 (0.10 ; 0.60)          2.00 (1.40 ; 2.50)
      20              \t Missing                        0                           0
              All (n = 100)
      1                    
      2         5.80 ± 0.95
      3  5.70 (5.00 ; 6.50)
      4  5.70 (4.30 ; 7.90)
      5                   0
      6                    
      7         3.20 ± 0.42
      8  3.20 (3.00 ; 3.42)
      9  3.20 (2.20 ; 4.40)
      10                  0
      11                   
      12        3.51 ± 2.10
      13 3.20 (1.50 ; 5.53)
      14 3.20 (1.00 ; 6.90)
      15                  0
      16                   
      17        1.14 ± 0.92
      18 1.00 (0.20 ; 2.00)
      19 1.00 (0.10 ; 2.50)
      20                  0

# Joining tables from nt_descriptive and nt_compare_tg

    Code
      as.data.frame(tab)
    Output
                        Variable Species: setosa (n = 50) Species: virginica (n = 50)
      1             Sepal.Length                                                     
      2             \t Mean ± SD              5.01 ± 0.35                 6.59 ± 0.64
      3  \t Median (Q25% ; Q75%)       5.00 (4.80 ; 5.20)          6.50 (6.23 ; 6.90)
      4    \t Median (Min ; Max)       5.00 (4.30 ; 5.80)          6.50 (4.90 ; 7.90)
      5               \t Missing                        0                           0
      6              Sepal.Width                                                     
      7             \t Mean ± SD              3.43 ± 0.38                 2.97 ± 0.32
      8  \t Median (Q25% ; Q75%)       3.40 (3.20 ; 3.68)          3.00 (2.80 ; 3.18)
      9    \t Median (Min ; Max)       3.40 (2.30 ; 4.40)          3.00 (2.20 ; 3.80)
      10              \t Missing                        0                           0
      11            Petal.Length                                                     
      12            \t Mean ± SD              1.46 ± 0.17                 5.55 ± 0.55
      13 \t Median (Q25% ; Q75%)       1.50 (1.40 ; 1.58)          5.55 (5.10 ; 5.88)
      14   \t Median (Min ; Max)       1.50 (1.00 ; 1.90)          5.55 (4.50 ; 6.90)
      15              \t Missing                        0                           0
      16             Petal.Width                                                     
      17            \t Mean ± SD              0.25 ± 0.11                 2.03 ± 0.27
      18 \t Median (Q25% ; Q75%)       0.20 (0.20 ; 0.30)          2.00 (1.80 ; 2.30)
      19   \t Median (Min ; Max)       0.20 (0.10 ; 0.60)          2.00 (1.40 ; 2.50)
      20              \t Missing                        0                           0
         p value
      1  < 0.001
      2         
      3         
      4         
      5         
      6  < 0.001
      7         
      8         
      9         
      10        
      11 < 0.001
      12        
      13        
      14        
      15        
      16 < 0.001
      17        
      18        
      19        
      20        

# Joining tables from nt_descriptive and nt_compare_mg

    Code
      as.data.frame(tab)
    Output
                        Variable Species: setosa (n = 50)
      1             Sepal.Length                         
      2             \t Mean ± SD              5.01 ± 0.35
      3  \t Median (Q25% ; Q75%)       5.00 (4.80 ; 5.20)
      4    \t Median (Min ; Max)       5.00 (4.30 ; 5.80)
      5               \t Missing                        0
      6              Sepal.Width                         
      7             \t Mean ± SD              3.43 ± 0.38
      8  \t Median (Q25% ; Q75%)       3.40 (3.20 ; 3.68)
      9    \t Median (Min ; Max)       3.40 (2.30 ; 4.40)
      10              \t Missing                        0
      11            Petal.Length                         
      12            \t Mean ± SD              1.46 ± 0.17
      13 \t Median (Q25% ; Q75%)       1.50 (1.40 ; 1.58)
      14   \t Median (Min ; Max)       1.50 (1.00 ; 1.90)
      15              \t Missing                        0
      16             Petal.Width                         
      17            \t Mean ± SD              0.25 ± 0.11
      18 \t Median (Q25% ; Q75%)       0.20 (0.20 ; 0.30)
      19   \t Median (Min ; Max)       0.20 (0.10 ; 0.60)
      20              \t Missing                        0
         Species: versicolor (n = 50) Species: virginica (n = 50) p value
      1                                                           < 0.001
      2                   5.94 ± 0.52                 6.59 ± 0.64        
      3            5.90 (5.60 ; 6.30)          6.50 (6.23 ; 6.90)        
      4            5.90 (4.90 ; 7.00)          6.50 (4.90 ; 7.90)        
      5                             0                           0        
      6                                                           < 0.001
      7                   2.77 ± 0.31                 2.97 ± 0.32        
      8            2.80 (2.52 ; 3.00)          3.00 (2.80 ; 3.18)        
      9            2.80 (2.00 ; 3.40)          3.00 (2.20 ; 3.80)        
      10                            0                           0        
      11                                                          < 0.001
      12                  4.26 ± 0.47                 5.55 ± 0.55        
      13           4.35 (4.00 ; 4.60)          5.55 (5.10 ; 5.88)        
      14           4.35 (3.00 ; 5.10)          5.55 (4.50 ; 6.90)        
      15                            0                           0        
      16                                                          < 0.001
      17                  1.33 ± 0.20                 2.03 ± 0.27        
      18           1.30 (1.20 ; 1.50)          2.00 (1.80 ; 2.30)        
      19           1.30 (1.00 ; 1.80)          2.00 (1.40 ; 2.50)        
      20                            0                           0        

