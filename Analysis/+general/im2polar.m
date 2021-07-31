% Convert an image to polar coordinates
%
% This function converts an image to polar coordinates. The function is useful to represent
% circular objects in image as though they are rectangular. One example is a transformation
% of a firing rate map recorded in a circular environment. Such image contains a circle of firing
% rates. Polar image will have degrees along y-axis and distance from the middle point of arena
% towards it's end along y-axis. If you map is [20 20]. Middle point is [10 10]. A pixel with
% coordinates [10 15] will have polar coordinates [0 5], i.e. 0 degrees, 5 is the distance from [10 10].
% There is no exact transformation and some values are linearly interpolated.
% The function creates an image that have 360 values along y-axis and appr. half of the width/height of im
% along the x-axis. If im is [36 36], then polar image size is [18 360].
%
%  USAGE
%    pim = general.im2polar(im)
%    im     2D matrix with the image to convert, usually a square matrix of size [NxN].
%    pim    Polar version of im. 2D matrix of size [N/2 360].
%
function pim = im2polar(im)
    [r, c] = size(im);
    cx = c/2 + 0.5;
    cy = r/2 + 0.5;

    rmax = min([cx-1, c-cx, cy-1, r-cy]);
    [theta, radius] = meshgrid(1:360, 0:rmax);

    xi = radius .* cosd(theta) + cx;  % Locations in image to interpolate data
    yi = radius .* sind(theta) + cy;  % from.

    [x, y] = meshgrid(1:c, 1:r);
    pim = interp2(x, y, im, xi, yi);
end