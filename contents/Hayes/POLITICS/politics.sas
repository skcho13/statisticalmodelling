data politics;
input pknow age educ sex income polint party libcon pdiscuss natnews npnews locnews talkrad ses news demoneg;
datalines;
    8   35   13     0   42.5     2    1    3    4    3    1  3.5  4.5  -.58 2.50 3.14
   12   40   14     0   42.5     2    2    6    3    0    4  3.5  3.5  -.34 2.50 2.85
   11   43   14     1  110.0     2    2    4    3    1    2  3.5  3.5   .50 2.16 2.00
   13   26   16     0  100.0     3    1    2    2    3    5  5.0  1.0   .85 4.33 2.00
   14   41   12     0   57.5     2    2    3    7    0    0  5.0  1.0  -.63 1.66 2.00
    8   41   12     1   80.0     2    1    3    3    0    7  3.5  1.0  -.34 3.50 2.00
   10   18   12     0   42.5     2    2    2    2    2    4  1.5  1.0  -.81 2.50 2.00
   15   31   14     1   90.0     3    2    7    7    2    0  2.0  4.5   .25 1.33 2.57
    9   18   11     0   30.0     3    1    5    5    3    7  7.0  1.0 -1.21 5.66 1.85
    6   72   12     1   30.0     2    3    6    7    0    7  7.0  1.0  -.97 4.66 2.57
    8   43   10     0   80.0     2    2    5    2    3    3   .5  1.0  -.82 2.16 2.57
    9   43   17     1   70.0     3    1    3    0    4    1  4.0  1.0   .71 3.00 1.85
   10   63   12     1   57.5     3    3    3    7    7    4  3.5  2.5  -.63 4.83 1.71
    2   32   14     0   57.5     2    2    5    7    3    1  1.5  1.0  -.15 1.83 2.42
    9   30   16     0   57.5     3    1    3    7    1    1  2.5  1.0   .32 1.50 1.71
    6   28   12     0   20.0     2    1    1    3    0    3  5.0  1.0 -1.10 2.66 2.00
    4   30   17     0   70.0     1    2    6    7    3    1  2.5  1.0   .71 2.16 3.14
   15   47   17     1  180.0     3    2    6    7    0    7   .0  3.0  2.09 2.33 2.42
    7   45   12     0   42.5     2    1    5    7    2    2  2.0  1.0  -.81 2.00 1.85
   12   60   14     0   42.5     2    1    3    5    0    7   .0  1.0  -.34 2.33 1.42
    6   52   14     0   42.5     3    3    5    7    5    3  7.0  3.5  -.34 5.00 1.57
    7   58   12     0   57.5     4    3    7    7    5    0  7.0  3.0  -.63 4.00 3.42
    9   38   16     0   80.0     3    2    5    4    2    0  2.0  1.0   .60 1.33 1.57
   11   58   13     0   42.5     2    3    3    3    0    7  3.5  1.0  -.58 3.50 1.71
    2   37   12     0   80.0     1    3    3    2    1    0  3.0  1.0  -.34 1.33 2.00
   20   54   14     1   42.5     4    2    7    4    7    0  1.5  3.0  -.34 2.83 2.57
    9   50   16     0   57.5     4    1    1    7    5    2  2.5  3.5   .32 3.16 3.00
   13   47   13     1   20.0     2    1    3    3    3    4  3.5  1.0  -.86 3.50 1.85
    7   26   14     0   42.5     2    1    2    5    2    2  3.5  1.0  -.34 2.50 2.57
    7   40   14     0   42.5     3    2    5    0    0    7  1.5  1.0  -.34 2.83 1.71
   12   24   15     0   30.0     2    2    5    7    2    0  4.5  2.0  -.26 2.16 3.20
   10   37   14     0  100.0     1    2    7    7    2    4  3.0  2.5   .37 3.00 2.71
   10   34   14     0   42.5     1    1    5    2    1    0  3.0  1.0  -.34 1.33 1.71
   14   56   16     1   70.0     3    2    6    7    7    0  2.0  1.0   .47 3.00 3.00
   12   25   16     1  120.0     2    1    3    5    4    3  3.0  1.0  1.10 3.33 1.85
   11   47   14     0   57.5     2    2    7    5    7    2  4.5  1.0  -.15 4.50 2.71
    9   57   11     1   30.0     3    1    3    6    3    4  3.5  3.5 -1.21 3.50 1.00
   20   49   17     1  120.0     4    2    6    7    0    7  3.5  4.0  1.34 3.50 2.42
    1   35   12     0   42.5     2    1    5    4    0    0   .5  1.0  -.81  .16 2.28
   10   51   13     0   57.5     2    3    5    7    7    7  7.0  1.0  -.39 7.00 1.33
   16   51   12     0  110.0     3    3    6    7    3    5  3.0  1.0   .02 3.66 2.85
   15   35   17     0   80.0     4    2    5    7    4    5  4.5  2.5   .84 4.50 3.00
   18   36   17     1   70.0     4    2    7    7    3    3  4.5  2.5   .71 3.50 2.71
   16   48   16     1   70.0     3    1    2    7    0    7   .0  3.0   .47 2.33 1.71
   16   40   16     0   42.5     4    2    5    7    0    2   .0  5.0   .13  .66 2.42
   13   33   15     1   57.5     3    1    5    3    5    2  3.0  3.5   .08 3.33 2.85
   16   29   17     0   30.0     4    1    2    5    4    0  1.0  1.0   .21 1.66 2.00
    2   20   12     0   42.5     2    2    4    3    4    1  1.0  1.0  -.81 2.00 2.85
   10   45   14     1  110.0     4    2    3    7    7    7  5.0  4.5   .50 6.33 3.00
   18   51   15     0  200.0     3    1    4    7    1    1   .5  1.0  1.86  .83 1.28
   16   65   14     1   57.5     3    2    6    7    7    2  7.0  3.0  -.15 5.33 2.00
   15   37   15     0   57.5     3    1    2    7    2    1  2.0  1.0   .08 1.66 1.14
   14   81   17     1   30.0     4    2    4    7    0    7  7.0  1.0   .21 4.66 3.71
    6   35   14     0   80.0     3    2    6    7    3    0  3.0  1.0   .12 2.00 2.85
   19   55   16     1   57.5     4    1    3    7    0    7   .0  1.0   .32 2.33 1.42
    3   22   12     0   20.0     2    3    5    6    0    1   .0  1.0 -1.10  .33 1.42
    6   74   12     0   20.0     4    1    3    0    7    7   .0  1.0 -1.10 4.66 1.83
    5   21   16     0   12.5     1    3    5    0    7    7  5.5  1.0  -.24 6.50 2.14
    6   76   12     0   20.0     3    2    5    7    7    6  7.0  1.0 -1.10 6.66 2.33
    4   50   17     0  150.0     2    2    5    5    1    2  2.5  1.0  1.71 1.83 2.42
   16   51   17     0  100.0     4    1    3    7    5    7   .5  1.0  1.09 4.16 1.28
    9   35   11     1   20.0     3    2    6    4    7    0  7.0  1.0 -1.34 4.66 2.60
    8   36   12     0   42.5     2    2    5    7    1    3   .0  1.0  -.81 1.33 2.28
    4   27   15     1   42.5     2    3    5    3    3    2   .0  1.0  -.10 1.66 3.42
   13   32   12     1   57.5     2    2    5    5    2    2  3.0  1.0  -.63 2.33 2.71
   10   67   14     1   42.5     2    2    4    0    1    7  5.0  1.0  -.34 4.33 2.33
    5   56   12     0   42.5     3    1    5    5    0    4  4.5  1.0  -.81 2.83 1.14
    5   46   16     0  100.0     3    2    2    7    7    1  5.0  1.0   .85 4.33 2.85
   15   63   16     1   57.5     4    3    6    7    7    7  6.0  5.0   .32 6.66 2.42
   14   46   12     0   90.0     2    1    2    3    3    0  3.5  1.0  -.22 2.16 1.28
   16   40   12     1   42.5     3    2    6    7    4    0  2.0  3.5  -.81 2.00 3.57
   18   41   16     0   80.0     4    3    5    7    7    0  4.0  4.5   .60 3.66 2.42
    7   32   15     1   70.0     3    3    2    7    2    1  2.5  3.0   .23 1.83 2.33
   17   54   16     1  100.0     4    2    6    7    7    7  3.5  2.5   .85 5.83 2.85
    9   30   14     1   30.0     3    2    2    0    0    4   .5  1.0  -.50 1.50 3.00
    9   62    8     0   57.5     3    1    2    0    7    2  2.0  4.0 -1.58 3.66 1.71
   10   53   16     1   57.5     3    1    5    7    3    3  4.0  1.0   .32 3.33 1.71
   13   56   13     0   42.5     4    2    6    7    7    7  7.0  4.5  -.58 7.00 2.00
   10   52   16     0   20.0     2    1    4    7    1    1   .5  1.0  -.14  .83 2.00
    7   37   13     0   80.0     2    1    5    7    0    2   .0  1.0  -.11  .66 3.00
   17   32   17     0    7.5     3    3    4    5    0    7   .5  1.0  -.06 2.50 2.71
   20   70   17     1   42.5     3    1    2    3    7    7  3.5  1.0   .37 5.83 1.14
   18   46   16     1   42.5     4    1    3    0    3    7  1.5  1.0   .13 3.83 1.83
    7   30   12     0   57.5     3    2    3    3    3    2  3.0  1.0  -.63 2.66 2.33
    8   26   15     1   30.0     3    2    4    3    1    1  7.0  3.5  -.26 3.00 2.14
    8   28   12     1   42.5     3    1    5    2    4    0  7.0  1.0  -.81 3.66 2.00
   11   65   16     0   42.5     3    1    7    7    0    7  7.0  1.0   .13 4.66 1.00
   17   35   16     1   57.5     4    2    5    7    3    7  3.5  4.5   .32 4.50 2.00
   12   54   13     1  110.0     2    1    3    0    0    0   .0  1.0   .26  .00 2.00
   15   28   17     1   57.5     2    2    4    7    0    7  2.0  1.0   .55 3.00 2.14
    5   82   12     0   20.0     4    1    3    0    7    1  7.0  3.5 -1.10 5.00 2.16
    7   90   13     0   20.0     4    3    6    0    7    0  2.5  3.5  -.86 3.16 2.00
   11   50   17     0   57.5     3    1    2    4    7    7  6.0  4.0   .55 6.66 1.85
    8   42   11     0    7.5     3    1    5    7    7    0  1.5  1.0 -1.49 2.83 1.33
   15   41   15     0   70.0     3    2    6    6    1    3  4.0  2.5   .23 2.66 3.00
   11   29   16     0   80.0     1    2    7    2    0    0   .5  3.0   .60  .16 2.85
   11   60   11     0   30.0     2    1    4    7    5    1  3.0  1.0 -1.21 3.00 1.57
    9   55   12     1   12.5     3    3    6    7    3    0  7.0  2.5 -1.19 3.33 1.71
   13   79   12     1   20.0     3    3    7    2    7    7  3.5  3.5 -1.10 5.83 2.57
    9   24   12     0   42.5     3    1    3    2    0    0  1.5  3.0  -.81  .50 3.42
   18   56   16     1  140.0     3    2    6    7    3    7  4.0  2.5  1.35 4.66 3.28
   14   64   17     1   30.0     4    3    7    7    7    7  6.0  1.0   .21 6.66 3.14
    7   50   12     0   20.0     2    3    5    7    0    0  3.5  1.0 -1.10 1.16 3.71
    6   42   11     1   20.0     3    3    5    6    7    2  3.5  1.0 -1.34 4.16 2.14
   15   78   12     1   30.0     3    2    6    0    4    7  2.5  2.5  -.97 4.50 2.28
    4   40   14     0   12.5     3    2    3    7    7    3  7.0  4.5  -.71 5.66 3.71
   18   27   17     0  110.0     3    1    1    7    0    5   .0  2.0  1.21 1.66 2.28
   14   41   17     0   42.5     3    2    6    7    2    7  3.0  3.5   .37 4.00 1.85
   15   74   15     0   30.0     3    1    3    7    7    7  5.5  1.0  -.26 6.50 1.85
   11   24   13     0   42.5     2    2    7    3    7    2  3.5  1.0  -.58 4.16 2.28
   15   43   17     1   70.0     4    1    1    7    7    7  1.0  1.0   .71 5.00 1.71
   15   22   16     1   12.5     2    1    2    7    2    4   .5  1.0  -.24 2.16 1.42
    8   44   12     0   57.5     2    2    5    7    2    1  5.5  1.0  -.63 2.83 2.85
    9   36   12     1   30.0     2    2    1    7    1    1  2.0  3.5  -.97 1.33 2.57
   13   34   17     0   57.5     3    2    6    2    0    0  1.5  1.0   .55  .50 2.00
   14   51   13     0   42.5     3    2    5    7    2    2  6.5  1.0  -.58 3.50 3.14
    5   35   13     0   30.0     2    2    6    3    0    0   .0  1.0  -.73  .00 3.33
   19   57   14     1   30.0     4    2    6    7    3    7  3.5  3.5  -.50 4.50 3.00
   14   75   14     0   70.0     3    1    5    7    7    7  7.0  1.0   .00 7.00 1.00
   11   24   11     1   30.0     3    1    5    2    1    7  7.0  1.0 -1.21 5.00 2.42
    8   64   15     1   12.5     4    1    6    7    7    7  7.0  1.0  -.48 7.00 1.00
   10   38   15     0   70.0     4    1    1    7    7    3  5.5  1.0   .23 5.16 1.71
    6   43   12     0   57.5     4    3    6    5    0    4   .0  1.0  -.63 1.33 2.71
   21   55   17     1   90.0     4    1    5    7    7    7  7.0  1.0   .96 7.00 2.00
    2   59    7     0   57.5     1    3    7    4    7    0  7.0  1.0 -1.82 4.66 1.42
    6   43   14     0   57.5     1    2    5    1    7    7  7.0  1.0  -.15 7.00 2.00
   14   74   16     1   42.5     4    1    3    7    4    5  1.5  3.5   .13 3.50 1.85
    6   40   17     1  100.0     2    2    6    2    1    7   .0  1.0  1.09 2.66 2.33
    1   39   12     1   42.5     2    3    5    0    0    0  1.0  1.0  -.81  .33 2.00
   14   39   17     1  100.0     3    2    6    7    2    7  2.5  1.0  1.09 3.83 2.85
   14   63   16     0   90.0     4    2    7    7    7    5  1.5  4.5   .72 4.50 3.14
   10   31   14     1   70.0     3    1    5    5    3    7  7.0  1.0   .00 5.66 1.85
   10   46   14     0   30.0     1    1    3    7    3    7  3.5  1.0  -.50 4.50 1.57
   10   72   16     0   57.5     2    2    7    5    5    6  2.5  1.0   .32 4.50 2.85
   11   70   12     0   42.5     3    2    5    7    7    1  3.5  1.0  -.81 3.83 2.33
   18   51   17     0   80.0     4    2    6    7    7    0  1.0  5.0   .84 2.66 2.71
   10   58   13     0   70.0     3    1    4    3    3    7  1.0  3.0  -.23 3.66 2.42
    8   37   14     0   42.5     2    1    3    7    1    1  4.0  2.5  -.34 2.00 2.85
   14   54   16     1   80.0     2    2    6    7    0    1   .0  4.0   .60  .33 2.66
   13   49   16     1   70.0     4    3    7    7    2    3  2.0  4.5   .47 2.33 4.00
    4   32   14     0  110.0     2    3    3    0    0    0  1.5  1.0   .50  .50 3.00
   19   49   16     1   80.0     3    1    3    4    0    7   .0  1.0   .60 2.33 1.42
   12   51   16     1    2.5     4    1    2    7    4    7  4.5  3.0  -.36 5.16 1.85
   13   24   15     0   70.0     4    1    3    7    7    7  3.5  3.0   .23 5.83 1.85
    4   21   11     0   42.5     2    2    3    0    7    7  3.5  1.0 -1.05 5.83 2.85
    9   46   14     1   42.5     4    1    5    7    4    5  4.5  1.0  -.34 4.50 1.85
    8   24   12     0   57.5     2    2    3    5    1    1  3.5  1.0  -.63 1.83 2.42
   12   32   16     0   80.0     2    2    7    7    2    0  5.0  4.0   .60 2.33 2.50
   12   46   12     0   70.0     3    2    7    7    7    6  3.0  1.0  -.47 5.33 3.00
   19   53   17     1  150.0     4    1    4    7    4    7  4.0  3.5  1.71 5.00 1.57
   11   59   15     0  140.0     4    3    5    7    7    1  3.5  1.0  1.11 3.83 2.33
   10   48   17     0   42.5     3    1    3    5    4    3  4.5  1.0   .37 3.83 1.42
   12   29   16     0   57.5     3    1    4    4    0    2   .5  2.5   .32  .83 1.83
   17   56   17     1   70.0     3    1    3    7    7    7  3.5  1.0   .71 5.83 2.00
    9   65   17     1   42.5     3    3    6    7    4    7  3.0  2.0   .37 4.66 2.66
   17   35   14     0   20.0     3    2    7    7    2    1  2.5  2.5  -.62 1.83 2.85
   14   45   12     1  120.0     4    2    7    7    0    0   .0  5.0   .15  .00 4.00
    9   53   16     0  100.0     3    2    6    7    5    7  7.0  3.0   .85 6.33 1.57
    7   56   12     0   57.5     1    1    5    0    0    7  4.0  1.0  -.63 3.66 2.00
    8   32   11     0   20.0     3    3    6    3    3    0  7.0  4.0 -1.34 3.33 2.85
    9   54   13     1   57.5     2    1    3    2    3    1  1.0  2.5  -.39 1.66 2.00
   16   49   12     1   70.0     3    3    3    7    2    2  7.0  4.5  -.47 3.66 2.14
   15   55   16     1   42.5     4    2    6    7    7    5  7.0  4.5   .13 6.33 2.00
   14   23   16     0  100.0     2    3    3    7    1    0  1.0  1.0   .85  .66 1.57
   12   59   12     0   80.0     3    1    2    7    3    2  1.5  1.0  -.34 2.16 1.00
   14   75   12     0   20.0     4    1    1    7    5    7  2.0  3.5 -1.10 4.66 1.42
    9   55   12     0   42.5     1    1    5    4    5    1  5.0  1.0  -.81 3.66 1.28
    8   43   14     0   57.5     1    1    1    7    0    2   .0  1.0  -.15  .66 2.50
   17   38   16     1  197.5     3    2    5    5    7    7   .0  1.0  2.07 4.66 2.00
   14   82   12     1   30.0     4    2    7    7    7    7  3.5  1.0  -.97 5.83 2.57
   13   37   16     1   70.0     3    2    5    7    1    2   .0  3.0   .47 1.00 1.85
   12   83   16     1   42.5     4    2    5    5    7    7  2.5  1.0   .13 5.50 2.28
   12   21   14     0   80.0     3    1    3    7    4    3   .5  2.5   .12 2.50 2.28
    7   19   12     0   20.0     2    1    4    3    0    1  1.5  1.0 -1.10  .83 1.42
   15   25   14     1   57.5     3    2    7    7    5    6  3.0  3.5  -.15 4.66 3.14
    7   41   12     1   80.0     1    2    2    7    0    2   .0  1.0  -.34  .66 2.66
   16   38   14     1   57.5     3    2    7    7    2    3  5.5  3.0  -.15 3.50 3.00
   14   35   17     1  120.0     4    2    5    7    7    7  7.0  3.5  1.34 7.00 3.42
   10   44   12     0   80.0     2    1    6    7    7    7  3.5  4.0  -.34 5.83 3.00
   11   35   14     1   80.0     2    2    2    2    3    0   .5  1.0   .12 1.16 2.28
   12   42   12     0  100.0     3    2    3    7    0    0  1.5  1.0  -.09  .50 3.60
    4   47   13     0   80.0     2    2    5    0    2    5   .0  1.0  -.11 2.33 1.66
    8   66   12     0   30.0     2    2    6    4    4    4  4.0  2.5  -.97 4.00 2.71
    6   30   16     0   30.0     3    3    3    7    1    5   .0  2.5  -.02 2.00 2.85
   19   46   17     1  150.0     4    1    3    7    5    4  4.5  1.0  1.71 4.50 1.42
    8   47   16     0   12.5     2    3    1    0    4    0  1.0  1.0  -.24 1.66 1.71
   12   60   16     0  130.0     3    2    6    7    5    7  2.5  1.0  1.22 4.83 2.57
    7   54   12     0   20.0     2    2    5    2    7    0  3.5  1.0 -1.10 3.50 3.00
   11   42   17     0  110.0     3    2    5    0    2    0  1.0  1.0  1.21 1.00 2.00
   19   49   16     1   70.0     4    1    3    7    0    7  1.0  1.0   .47 2.66 2.00
   12   38   12     0   57.5     2    2    5    1    1    1   .5  4.0  -.63  .83 2.14
    7   36   12     1   42.5     3    1    5    2    3    1  2.0  1.0  -.81 2.00 1.71
    9   52   13     1   70.0     3    1    5    0    7    3  7.0  1.0  -.23 5.66 1.71
   15   44   14     1   70.0     3    1    4    4    4    7  2.0  3.5   .00 4.33 2.00
   12   35   16     1   80.0     3    1    5    7    4    7  5.0  4.0   .60 5.33 1.00
   19   53   17     1   57.5     4    1    2    7    4    7  5.0  1.0   .55 5.33 1.57
    8   68   12     1   20.0     3    1    5    3    2    7  2.5  1.0 -1.10 3.83 2.00
   14   46   14     0   42.5     3    2    7    3    3    6   .0  1.0  -.34 3.00 2.57
   14   43   11     1   57.5     1    1    5    5    1    2  1.5  1.0  -.87 1.50 1.14
   11   54   12     1   90.0     4    1    6    7    5    5  5.0  1.0  -.22 5.00 2.28
   17   38   16     1  140.0     4    2    5    7    7    7  2.0  2.5  1.35 5.33 3.42
   14   84   17     1  100.0     3    1    3    0    2    3  4.0  1.0  1.09 3.00 1.85
    6   64   16     0    2.5     3    1    6    0    7    3  5.0  1.0  -.36 5.00 1.83
   10   57   14     0   57.5     3    1    5    7    2    2  1.0  1.0  -.15 1.66 2.42
    9   47   16     0   30.0     3    3    1    7    7    7  7.0  1.0  -.02 7.00 1.71
   14   29   17     1   12.5     3    2    6    3    0    0  3.0  1.0   .00 1.00 2.57
   15   48   16     1   70.0     2    1    5    7    7    7  3.5  1.0   .47 5.83 1.71
   12   58   16     1  130.0     2    2    6    7    3    7  3.0  1.0  1.22 4.33 1.85
   20   47   15     1   80.0     4    1    3    5    5    7  4.0  3.5   .36 5.33 1.83
    7   51   12     0   70.0     1    1    4    4    7    7  7.0  1.0  -.47 7.00 1.71
    7   36   14     0   80.0     2    2    3    0    0    0  7.0  3.0   .12 2.33 2.85
   11   19   11     1   42.5     3    1    2    4    1    0  1.0  1.0 -1.05  .66 2.14
   15   52   17     1  180.0     2    2    6    2    7    1  1.5  3.0  2.09 3.16 2.71
    8   73   12     0   20.0     3    2    5    0    7    7  4.0  1.0 -1.10 6.00 1.50
    3   36   12     0   57.5     1    1    6    0    0    0  1.0  1.0  -.63  .33 2.14
   16   35   17     0  180.0     3    2    5    7    7    1  2.0  3.0  2.09 3.33 2.14
   18   27   17     1   57.5     4    2    6    4    0    7  1.0  2.0   .55 2.66 2.85
   16   29   17     0  130.0     3    2    6    3    2    0   .0  3.0  1.46  .66 2.85
   14   72   14     0   30.0     3    2    7    7    6    5  6.0  2.0  -.50 5.66 3.00
   14   45   13     1  120.0     4    1    3    7    7    7  3.5  3.0   .38 5.83 2.85
    8   38   11     1   42.5     1    1    1    7    7    0  3.5  3.5 -1.05 3.50 2.00
   12   27   17     1   57.5     3    1    4    2    0    2   .0  3.5   .55  .66 2.00
   10   58   13     0   57.5     3    2    7    7    0    5   .0  1.0  -.39 1.66 3.42
   18   31   17     1  130.0     4    2    5    3    7    7  1.5  2.0  1.46 5.16 2.42
    5   37   12     1   20.0     3    3    5    4    2    1  2.0  1.0 -1.10 1.66 1.85
    8   42   16     0   57.5     3    2    5    1    1    3  1.0  3.5   .32 1.66 2.00
    8   52   12     1   42.5     3    1    5    3    7    0  7.0  1.0  -.81 4.66 1.66
   15   61   16     0  150.0     2    2    5    7    2    0  1.0  4.5  1.48 1.00 2.42
   12   30   17     0   30.0     3    2    6    7    0    2  1.0  1.0   .21 1.00 2.14
   19   48   17     1   80.0     4    1    1    7    0    7  1.0  1.0   .84 2.66 1.14
   19   33   17     1  180.0     4    3    5    4    7    5  3.5  2.5  2.09 5.16 2.14
   12   42   16     1   30.0     4    3    5    2    1    3  1.0  4.0  -.02 1.66 1.71
   11   57   14     0   70.0     4    2    5    7    7    7  5.0  1.0   .00 6.33 2.14
   14   24   14     0   12.5     2    2    5    7    3    4  1.5  1.0  -.71 2.83 2.57
   13   49   17     0   42.5     4    1    2    7    5    3  4.0  3.5   .37 4.00 1.28
    6   56   15     1   30.0     4    1    5    5    3    2  2.5  2.0  -.26 2.50 1.85
    2   28   13     1   30.0     1    3    5    2    1    2  2.5  1.0  -.73 1.83 2.42
   13   46   16     0   90.0     2    2    5    7    3    7  1.5  2.0   .72 3.83 2.28
   14   33   14     1   57.5     3    2    6    7    3    1  2.0  3.0  -.15 2.00 3.00
   18   45   16     1   80.0     4    2    5    7    7    7  5.0  1.0   .60 6.33 2.14
   11   25   12     0   57.5     1    1    3    0    2    1  4.0  1.0  -.63 2.33 1.85
    9   53   13     0   42.5     3    3    3    0    7    1  2.5  5.0  -.58 3.50 1.71
   11   56   12     1   30.0     3    2    7    7    3    7   .0  3.5  -.97 3.33 4.00
    9   56   14     1   42.5     3    1    3    7    7    1  5.5  3.0  -.34 4.50 1.16
   14   69   12     0   20.0     4    1    2    0    7    0  7.0  4.5 -1.10 4.66 1.42
    3   21   15     0   70.0     2    1    2    1    0    1   .0  1.0   .23  .33 2.00
   15   63   12     1   42.5     3    2    6    4    7    7  3.5  3.5  -.81 5.83 3.71
    6   38   15     0    2.5     3    2    5    7    0    0   .0  1.0  -.60  .00 2.85
   16   50   12     0  200.0     3    2    6    7    2    1   .0  4.0  1.15 1.00 2.57
   13   52   14     0   70.0     3    2    7    5    3    0  1.5  4.0   .00 1.50 2.71
    8   28   13     0   30.0     1    1    3    2    0    2  3.5  1.0  -.73 1.83 1.71
   13   33   15     1   90.0     3    2    6    7    0    7  1.5  2.5   .49 2.83 2.85
    9   58   13     1   42.5     1    3    5    0    7    7  7.0  1.0  -.58 7.00 2.14
   17   32   16     1  100.0     4    2    6    7    3    0  2.0  1.0   .85 1.66 2.71
   14   40   17     0   42.5     3    1    2    7    5    2  5.5  1.0   .37 4.16 1.71
   10   41   14     0   57.5     3    2    6    7    7    3  2.0  4.0  -.15 4.00 2.85
   20   32   16     0   57.5     3    1    1    7    2    4   .0  4.0   .32 2.00 1.42
    6   32   15     0   42.5     2    3    4    7    0    1  4.5  1.0  -.10 1.83 2.57
   15   58   13     1   42.5     4    2    7    2    0    7   .0  3.5  -.58 2.33 2.57
   18   33   17     1  110.0     4    2    4    7    7    5  6.0  4.5  1.21 6.00 2.57
   17   51   16     1   80.0     4    2    6    3    3    7   .0  1.0   .60 3.33 2.42
   17   32   16     1   57.5     4    2    7    7    2    2  1.5  3.5   .32 1.83 3.00
    0   24   16     1  100.0     2    2    7    7    0    3   .0  1.0   .85 1.00 2.85
    8   57   14     1   70.0     2    1    5    7    4    0  4.0  1.0   .00 2.66 2.00
   15   50   17     0   80.0     3    1    2    4    0    4   .0  1.0   .84 1.33 2.00
    8   46   14     0   57.5     3    3    3    3    7    5  1.0  1.0  -.15 4.33 2.20
    6   37   13     0   70.0     2    1    6    2    0    2  3.0  1.0  -.23 1.66 2.00
   12   28   14     1   30.0     3    2    6    4    5    2  3.5  3.0  -.50 3.50 1.71
    5   42   14     1   57.5     2    1    2    4    4    7  7.0  1.0  -.15 6.00 1.42
    8   40   15     1   80.0     2    1    2    4    6    4  5.5  1.0   .36 5.16 2.16
   13   54   16     0  160.0     3    2    7    7    7    7   .5  3.0  1.60 4.83 2.71
   16   39   17     1   30.0     4    3    2    1    0    7   .0  1.0   .21 2.33 2.00
   13   45   17     1  160.0     4    2    7    7    7    2   .0  3.5  1.84 3.00 3.85
   12   73   12     0   20.0     4    3    6    4    5    7  4.5  4.5 -1.10 5.50 2.50
   12   37   16     0  100.0     3    3    5    3    5    1  2.5  1.0   .85 2.83 2.71
    7   46   16     0   80.0     2    1    3    7    7    7  3.0  5.0   .60 5.66 1.71
   17   30   15     1   42.5     3    1    3    7    6    3   .0  1.0  -.10 3.00 1.71
   14   56   15     1   70.0     2    1    7    7    7    7  3.5  3.0   .23 5.83 1.85
   14   35   14     1  200.0     4    1    1    3    1    5   .0  1.0  1.63 2.00 2.14
    6   75   14     0   30.0     3    1    5    0    0    7  7.0  1.0  -.50 4.66 2.14
   15   62   16     1   57.5     4    2    6    7    5    7  3.5  4.5   .32 5.16 3.00
    8   41   12     0   30.0     3    1    4    7    0    3  3.5  1.0  -.97 2.16 1.83
   16   55   17     1  120.0     2    1    3    3    3    2  2.0  4.0  1.34 2.33 3.14
   10   39   14     0   30.0     2    1    3    0    3    2  3.5  1.0  -.50 2.83 2.42
   14   43   16     1  100.0     3    1    3    4    4    1  2.5  1.0   .85 2.50 1.80
    7   21   14     1   30.0     2    3    2    4    0    3  3.5  3.0  -.50 2.16 2.00
    9   51   16     0   70.0     4    1    6    7    7    7  7.0  3.5   .47 7.00 1.57
   15   40   12     1   80.0     3    1    6    3    4    7  1.5  3.0  -.34 4.16 1.57
   18   28   16     0   70.0     3    1    3    7    7    7  1.5  3.5   .47 5.16 2.42
   15   53   13     1  130.0     3    2    6    7    4    7  1.0  1.0   .51 4.00 2.14
    3   35   16     0   57.5     2    2    6    2    0    0   .5  1.0   .32  .16 2.00
   12   54   14     1   57.5     3    1    3    3    3    7  3.5  1.0  -.15 4.50 2.28
    8   55   16     1   57.5     2    1    5    7    5    0  3.0  2.0   .32 2.66 2.57
   11   23   12     1   42.5     2    1    5    0    4    7  7.0  3.5  -.81 6.00 2.00
    9   30   16     1   57.5     2    3    5    3    3    2   .0  4.0   .32 1.66 1.85
   14   42   16     1   80.0     3    1    3    7    4    5  2.0  3.0   .60 3.66 2.14
    2   54    7     1   12.5     2    1    5    0    0    0  2.5  1.0 -2.38  .83 3.40
    7   43    6     1   12.5     1    3    5    0    2    0  3.5  1.0 -2.62 1.83 2.16
   11   40   16     1   90.0     3    1    4    0    3    3   .0  1.0   .72 2.00 1.00
   14   22   13     0   20.0     2    3    2    7    2    0  5.0  1.0  -.86 2.33 1.57
   11   50   13     1  110.0     3    2    5    7    3    3  2.0  1.0   .26 2.66 2.14
   11   36   13     1   80.0     3    2    6    7    7    1  3.0  4.0  -.11 3.66 3.57
   15   49   13     1  200.0     4    2    7    7    5    7  3.0  4.5  1.39 5.00 3.57
   11   37   12     0  110.0     2    2    2    7    1    7  1.0  1.0   .02 3.00 1.85
   11   27   15     0   42.5     1    3    6    2    1    0  5.0  1.0  -.10 2.00 2.00
   10   42   15     1  110.0     2    3    3    3    3    2  1.0  1.0   .74 2.00 1.83
   16   37   17     1   42.5     2    2    5    7    5    0  3.0  4.0   .37 2.66 1.85
    8   71    8     1   20.0     2    1    3    0    7    0  3.5  1.0 -2.05 3.50 2.85
   11   38   16     1  200.0     3    2    5    2    7    7   .0  1.0  2.10 4.66 1.71
   13   42   16     0  170.0     2    2    6    7    3    7  4.5  1.0  1.73 4.83 2.85
   12   33   16     0   42.5     4    1    2    5    3    0  1.5  1.0   .13 1.50 2.00
    8   28   14     1   42.5     2    2    3    3    3    0  4.5  1.0  -.34 2.50 2.42
    8   29   12     0   42.5     2    2    2    0    0    0  2.0  1.0  -.81  .66 2.28
   12   67   13     0   12.5     2    3    5    0    7    7  2.0  1.0  -.95 5.33 1.71
   10   33   16     0   70.0     4    2    5    7    3    2  2.0  4.0   .47 2.33 2.28
    6   36   12     0   30.0     3    1    5    3    3    0  3.0  4.0  -.97 2.00 2.14
    8   38   16     0   30.0     2    2    6    3    0    0   .5  1.0  -.02  .16 2.85
    8   41   13     0   57.5     3    2    7    7    0    3  1.0  3.0  -.39 1.33 2.28
   15   45   15     0   42.5     4    1    3    7    7    7  5.5  2.5  -.10 6.50 2.00
    6   39   12     1   30.0     2    1    7    4    3    2  7.0  3.0  -.97 4.00 1.57
   19   55   16     1   70.0     3    3    5    7    1    7  2.0  1.0   .47 3.33 1.42
    9   61   12     1   90.0     4    1    6    7    7    7  7.0  1.0  -.22 7.00 3.28
   18   28   16     1   42.5     4    1    3    7    4    4  3.0  3.5   .13 3.66 1.71
   17   53   16     0  170.0     4    1    3    4    4    3  1.5  1.0  1.73 2.83 1.71
   16   37   16     1  100.0     3    1    3    7    0    7  1.0  1.0   .85 2.66 2.00
   15   44   14     1   57.5     4    2    7    7    7    0  3.5  2.5  -.15 3.50 3.42
   13   68   13     1   80.0     3    2    6    2    0    2   .0  3.0  -.11  .66 2.42
   15   60   17     1   30.0     3    2    7    4    1    5  1.5  1.0   .21 2.50 2.42
   13   32   16     1  100.0     2    1    2    3    2    3  1.5  1.0   .85 2.16 2.00
    7   23   16     1   30.0     2    2    3    3    0    5  2.0  1.0  -.02 2.33 2.00
   11   75   14     0   30.0     4    2    7    7    7    7  7.0  2.5  -.50 7.00 2.71
   12   49   17     0   90.0     4    1    2    7    2    7  1.0  3.0   .96 3.33 1.80
   12   36   16     1   90.0     2    1    5    3    0    0  1.5  1.0   .72  .50 1.33
   15   56   17     0   57.5     4    2    5    7    7    7  3.5  3.5   .55 5.83 2.57
   20   28   16     0   70.0     3    1    3    7    7    0   .0  1.0   .47 2.33 1.28
   14   57   17     1  170.0     4    2    7    7    7    4  3.5  3.5  1.96 4.83 2.85
    5   51   12     0   57.5     2    2    5    1    0    0   .0  1.0  -.63  .00 2.83
    3   58   13     0   57.5     2    1    5    7    4    7  4.5  3.5  -.39 5.16 2.00
   15   23   13     1   90.0     3    2    5    7    3    0  2.0  4.0   .01 1.66 1.42
   10   73   12     0   12.5     2    1    5    6    7    7  7.0  1.0 -1.19 7.00 1.57
run;
