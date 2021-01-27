function [cmap]=makecolormap(namemap,nsteps)
%makecolormap this function returns a colormap
% adapted from ColorBrewer, and possibly interpolated.
% [cmap]=makecolormap(namemap,nsteps)
% returns a colormap in nsteps
% [map]=makecolormap(namemap) takes nsteps = 100
% namemap = 'redwhiteblue', 'redwhite', 'redwhitegreen',
% 'brownwhitegreenblue', 'brownwhitegreen','whitegreen','redwhiteblue2',
% 'redwhite', 'redwhite2','redyellowgreen', 'qualitative1', 'qualitative2',
% 'qualitative3'
% Qualitative colormaps are not interpolated
%
% written by Taco Broerse, 2020

narginchk(1,2);

if nargin==1
    nsteps=100;
end

if strcmp(namemap,'redwhiteblue')
    % colorbrewer default
    cmap=[103 0 31;
        178 24 43;
        214 96 77;
        244 165 130;
        253 219 199;
        247 247 247;
        209 229 240;
        146 197 222;
        67 147 195;
        33 102 172;
        5 48 97];
    
elseif strcmp(namemap,'redwhitegreen')
    cmap=[103 0 31;
        178 24 43;
        214 96 77;
        244 165 130;
        253 219 199;
        247 247 247;
        217   239   139;
        166   217   106;
        102   189    99;
        26   152    80;
        0   104    55];
    
    
elseif strcmp(namemap,'brownwhitegreenblue')
    cmap=[
        84    48     5
        140    81    10
        191   129    45
        223   194   125
        246   232   195
        245   245   245
        199   234   229
        128   205   193
        53   151   143
        1   102    94
        0    60    48];
elseif strcmp(namemap,'brownwhitegreen')
    cmap=[
        102    37     6
        153    52     4
        204    76     2
        236   112    20
        254   153    41
        254   196    79
        254   227   145
        255   247   188
        255   255   229
        247   252   253
        229   245   249
        204   236   230
        153   216   201
        102   194   164
        65   174   118
        35   139    69
        0   109    44
        0    68    27];
elseif strcmp(namemap,'brownwhitegreen2')
    cmap =[

    84    48     5
   140    81    10
   191   129    45
   223   194   125
   246   232   195
   245   245   245
   199   234   229
   128   205   193
    53   151   143
     1   102    94
     0    60    48];
elseif strcmp(namemap,'whitegreen')
    cmap=[
        
    0    68    27
    0   109    44
    35   139    69
    65   174   118
    102   194   164
    153   216   201
    204   236   230
    229   245   249
    247   252   253
    255   255   255];
elseif strcmp(namemap,'redwhiteblue2')
    % composed of red to white and white to blue from colorbrewer
    cmap=[
        103     0    13
        165    15    21
        203    24    29
        239    59    44
        251   106    74
        252   146   114
        252   187   161
        254   224   210
        255   255   255
        222   235   247
        198   219   239
        158   202   225
        107   174   214
        66   146   198
        33   113   181
        8    81   156
        8    48   107];
elseif strcmp(namemap,'redwhite')
    cmap=[103 0 31;
        178 24 43;
        214 96 77;
        244 165 130;
        253 219 199;
        247 247 247];
elseif strcmp(namemap,'redwhite2')
    cmap=[
        103     0    13
        165    15    21
        203    24    29
        239    59    44
        251   106    74
        252   146   114
        252   187   161
        254   224   210
        255   255   255];
elseif strcmp(namemap,'redyellowgreen')
    cmap=[215 48 39;
        244 109 67;
        253 174 97;
        254 224 139;
        255 255 191;
        217 239 139;
        166 217 106;
        102 189 99;
        26 152 80];
elseif strcmp(namemap,'qualitative1')
    nsteps=8;
    cmap=[27,158,119;
        217,95,2;
        117,112,179;
        231,41,138;
        102,166,30;
        230,171,2;
        166,118,29;
        102,102,102];
    cmap=flip(cmap);
elseif strcmp(namemap,'qualitative2')
    nsteps=12;
    cmap=[
        166   206   227;
        31   120   180;
        178   223   138;
        51   160    44;
        251   154   153;
        227    26    28;
        253   191   111;
        255   127     0;
        202   178   214;
        106    61   154;
        255   255   153;
        177    89    40];
    cmap=flip(cmap);
elseif strcmp(namemap,'qualitative3')
    nsteps=8;
    cmap=[228,26,28
        55,126,184
        77,175,74
        152,78,163
        255,127,0
        255,255,51
        166,86,40
        247,129,191];
    cmap=flip(cmap);
else
    error('unknown color map')
end
% normalise
cmap=(flip(cmap)/255);



%% interpolate colormap
cmap=interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,nsteps));

end