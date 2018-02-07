function [d_min, d_max] = max_min_diameter(BW)

imsize = size(BW);
BW(1:5, 1:5) = zeros(5, 5);

find_white = 0;
i = 1;
while ((find_white == 0) && (i <= imsize(1)))
    j = find(BW(i, :) == 1, 1, 'first');
    if (min(size(j)) ~= 0)
        find_white = 1;
        top_i = i;
        top_j = j;
    else
        i = i + 1;
    end;
end;

find_black = 0;
while ((find_black == 0) && (i <= imsize(1)))
    j = find(BW(i, :) == 1, 1, 'first');
    if (min(size(j)) == 0)
        find_black = 1;
        j = find(BW(i-1, :) == 1, 1, 'first');        
        bottom_i = i -1;
        bottom_j = j;
    else
        i = i + 1;
    end;
end;

if (i >= imsize(1))
    i = imsize(1);
    j = find(BW(i, :) == 1, 1, 'first');
    bottom_i = i;
    bottom_j = j;
end;

find_white = 0;
j = 1.0;
while ((find_white == 0) && (j <= imsize(2)))
    i = find(BW(:, j) == 1, 1, 'first');
    if (min(size(i)) ~= 0)
        find_white = 1;
        right_i = i;
        right_j = j;
    else
        j = j + 1;
    end;
end;

find_black = 0;
j = ceil( imsize(2) /2); 
while ((find_black == 0) && (j <= imsize(2)))
    i = find(BW(:, j) == 1, 1, 'first');
    if (min(size(i)) == 0)
        find_black = 1;
        i = find(BW(:, j-1) == 1, 1, 'last');        
        left_i = i;
        left_j = j -1;
    else
        j = j + 1;
    end;
end;

if (j >= imsize(2))
    j = imsize(2);
    i = find(BW(:, j) == 1, 1, 'last');
    left_i = i;
    left_j = j;
end;

d1 = norm([top_i - bottom_i, top_j - bottom_j]);
d2 = norm([right_i - left_i, right_j - left_j]);

%disp([top_i, top_j, bottom_i, bottom_j, right_i, right_j, left_i, left_j])

d_min = min([d2, d1]);
d_max = max([d2, d1]);
