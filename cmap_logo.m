function cmap_logo(h_axes)

[I, ~, alpha] = imread('cMapApp logo small.png');
alpha(all(I == 255, 3)) = 0;

h_logo = imshow(uint8(I), 'Parent', h_axes);

set(h_logo, 'AlphaData', alpha);
axis(h_axes, 'square');
