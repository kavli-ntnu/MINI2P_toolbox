function output=circularSmooth(input, smoothnumber)

input_extent=[input;input(1:1:smoothnumber)];
input_extent=smooth(input_extent,smoothnumber);
output=input_extent(1:size(input,1));
end