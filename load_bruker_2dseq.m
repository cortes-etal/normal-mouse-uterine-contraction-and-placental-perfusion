function [ image_space ] = load_bruker_2dseq( num_y, filename, data_type )

%*****************************************************************************
%
%  load_bruker_2dseq.m
%
%  Description:  Automatically or interactively identifies a data file, 
%    loads it, and outputs image_space data file to be further manipulated.
%
%  Output arguments:
%    image_space - (2-D) image data created by the ParaVision software
%
%  Input arguments:
%    num_y - (1x1) number of points along the y direction (int)
%
%    filename - (1xn) optional filename (str)
%
%    data_type - optional data type string (e.g., int32). Defaults to int16
%
%  Functions calling this function:
%    none
%
%  Functions called by this function:
%    none
%
%  Usage:
%    [ image_space ] = load_bruker_2dseq( num_y );
%
%    [ image_space ] = load_bruker_2dseq( num_y, filename );
%
%    [ image_space ] = load_bruker_2dseq( num_y, filename, data_type );
%
%*****************************************************************************

% --- Initialize variables.
image_space = [];

% --- Let the user choose a filename if one is not given.
if ( nargin < 2 )
  [ fn, pn ] = uigetfile( { '*', 'All Files' }, ...
    'Open a Bruker 2dseq file.' );
  if ( isequal( fn, 0 ) )
    return;
  end; % if
  data_type='int16';
elseif ( nargin < 3 )
  data_type='int16';
else
  [ pn, fn ] = fileparts( filename );
  if ( isempty( pn ) )
    pn='.';
  end % if
end % if-else

% --- Open the file.
fid_in = fopen( fullfile( pn, fn ), 'r', 'l' );

if ( fid_in == -1 )
  disp( [ 'ERROR: Cannot open file: ', filename ] );
  return;
end; % if

% --- Read the entire file into one vector.
file_stream = double( fread( fid_in, data_type ) );

% --- Close the file.
fclose( fid_in );

% --- Reformat the vector into a 2-D matrix of y-x.
image_space = reshape( file_stream, num_y, length( file_stream ) / num_y );

% >>> EOF

