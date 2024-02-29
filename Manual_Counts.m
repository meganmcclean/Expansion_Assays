%MANUAL_COUNTS
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax: Manual_Counts
%
% Inputs:
%   None
%   The user selects files to be analyzed manually in respose to a
%   uigetfile call. 
%   Expectations for the file format, naming, and organization:
%   TBD
%
% Outputs:
%    TBD
%
% Example: 
%   TBD
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TBD

% Author: FirsMegan N. McClean
% Department of Biomedical Engineering, University of Wisconsin Madison
% Engineering Centers Buidlgin 3156 1550 Engineering Drive
% Madison, WI 53705
% email: mmcclean@wisc.edu
% Website: https://github.com/meganmcclean/Expansion_Assays
% February 2024; Last revision: February 2024
%
% Acknowledgements: Most of this code was recycled from code written by
% Kevin Stindt, Ph.D. to automatically perform analysis of expansion assay
% colony sectors. 


close all; clear all; 

%% Select all files to load
% Make sure to select ALL files of a particular day/sample (i.e. GFP, YFP, and
% RFP image)
[file, folder] =  uigetfile('.tif','Select files','MultiSelect','on');


%% Convert file info into struct (fancy table)
data = [];
data0 = table();
for i = 1:length(file)
    C = strsplit(file{i},{'_','.'});
    data0.folder = cellstr(folder);
    data0.filename = file(i);
    data0.Cooperator = categorical(C(1));
    data0.Day = categorical(C(2));
    data0.Carbon = categorical(C(3));
    data0.Ratio = categorical(C(4));
    data0.Fluor = categorical(C(5));
    data0.Rep = categorical(C(6));
    data = vertcat(data,data0);
end  

%% Convert Days (could do this in Finder too...)
data.Day(data.Day=='D2') = '2';
data.Day(data.Day=='D3') = '3';
data.Day(data.Day=='D4') = '4';
data.Day(data.Day=='D5') = '5';
data.Day(data.Day=='D6') = '6';
data.Day(data.Day=='D7') = '7';
data.Day = double(string(data.Day));

data.Rep(data.Rep=='R1') = '1';
data.Rep(data.Rep=='R2') = '2';
data.Rep(data.Rep=='R3') = '3';
data.Rep = double(string(data.Rep));

data.Rat_Cheat(data.Ratio=='1-1') = 1;
data.Rat_Cheat(data.Ratio=='1-30') = 1;
data.Rat_Cheat(data.Ratio=='30-1') = 30;
data.Rat_Coop(data.Ratio=='1-1') = 1;
data.Rat_Coop(data.Ratio=='1-30') = 30;
data.Rat_Coop(data.Ratio=='30-1') = 1;

%% Iterate through the samples allowing the user to manually measure what they want to measure

% Setup radii measurements to use for all images to be consistent
colony_radius=160; %Radius picked by eye to make sure that all of the colony is included


    radii=input("Enter a vector of radial distances at which you want to measure your factor of interest (don't forget square bracket):");
%catch
    if isempty(radii) %If the user doesn't enter radii, they are stuck with these.
        radii=[60 80 100 110 120 130 140 150 160];
    end
%end

% Force the last point in radii to be colony_radius
if radii(end)~=colony_radius
    radii=[radii colony_radius];
end

data_sectors=[]; %To make a new table with only the sector information
data_clicks=[]; %To make a new table with the click information (including x,y positions)

mkdir Image_Output

for j = 1:size(data,1)
    
    % some cleanup
    clearvars -except data data_sectors data_clicks file folder j colony_radius radii
    close all; % close all figures

    d = data.Day(j);
    coop = data.Cooperator(j);
    carb = data.Carbon(j);
    rat = data.Ratio(j);
    rep = data.Rep(j);
     
    jR = find(data.Day==d & data.Cooperator==coop & data.Carbon==carb &...
            data.Ratio==rat & data.Rep==rep & data.Fluor=='RFP');
    jY = find(data.Day==d & data.Cooperator==coop & data.Carbon==carb &...
            data.Ratio==rat & data.Rep==rep & data.Fluor=='YFP');
    
    if data.Fluor(j)=="GFP" %Cycle through samples once, arbitrarily choosing GFP channel to make sure you don't double up

        %Import each Tiff file, binarize, fill, for both green then red channels
        pathG = [data.folder{j} data.filename{j}];
        pathR = [data.folder{jR} data.filename{jR}];
        pathY = [data.folder{jY} data.filename{jY}];
        
        % Load image and convert to greyscale
        im = imread(pathY);
        % Debugging
        %imshow(im,[])
        %pause

        % Find colony edge from binarized image
        im_binary = imadjust(im);
        im_binary = imbinarize(im_binary);
        % Used center of mass from regionprops to find center of the YFP
        % image (where all cells are labelled)
        stats=regionprops('table',im_binary, 'Centroid','MajorAxisLength','MinorAxisLength','Area');
        stats=sortrows(stats, 'Area','descend');
        colony_center=stats.Centroid(1,:);
       
         
        %Add relevant information to the data table
         data.Center{j} = colony_center;
         data.Center{jR}=colony_center;
         data.Center{jY}=colony_center;
         data.MaxRadiusPixels(j) = colony_radius;
         data.MaxRadiusPixels(jR) = colony_radius;
         data.MaxRadiusPixels(jY) = colony_radius;
  
        % Center the image on the central yeast spot
        spot_center=colony_center;
        spot_radius=colony_radius;
   
        image_center = [size(im,2)/2, size(im,1)/2];
        im = imtranslate(im,image_center - spot_center);
        colony_center = colony_center + (image_center - spot_center);
        
        % debugging
        % close all
        % imshow(im,[]); hold;

        %Read in the GFP image
        im2=imread(pathG);
        im2=imtranslate(im2, image_center-spot_center);
        
            
        % Display the GFP image (cheaters are bright, cooperators are dark)
        % Display circles (red) at the user picked radii
        % plot the center as a red dot
        imshow(im2,[]); hold; 
        viscircles(repmat(colony_center,[length(radii),1]), radii, 'Linewidth',0.25)
        plot(colony_center(1,1), colony_center(1,2), 'ro')
        
        
        % For each user selected radii, user will count just inside of that
        % radius
        for r_i=1:length(radii) 
            temp_table=data(j,:);
            temp_table_clicks=data(j,:);
            temp_table.radii=radii(r_i);
            temp_table_clicks.radii=radii(r_i);
            viscircles(colony_center,radii(r_i),'Color','blue', 'Linewidth',0.25)
            disp(strcat('Pick cheater (GFP) sector to start with. Will be marked by a *. Radius is R=',num2str(radii(r_i))))
            try  %Check skipping here
                gtext('*','FontSize',14)
                disp('Count cheaters (GFP) with left mouse click. Count cooperator (dark) with a right mouse click. NOTE: Click the cheater you picked to start with to START')
                [x,y, Button]=ginput();
                plot(x,y,'c.')
               
                %Keep clicked points and button information
                temp_table_clicks.x{1} = x;
                temp_table_clicks.y{1} = y;
                temp_table_clicks.Button{1} = Button;
               
                
                temp_table.Flag=0; %Defalt flag value is 0= no problems
                
                %Ask the user if they want a different flag
                try
                    temp_table.Flag=input("Points you counted are displayed. Enter a flag if needed (integer >0)");
                catch
                end   
                temp_table_clicks.Flag=temp_table.Flag;
             
                temp_table.left_click=length(find(Button==1));
                temp_table.right_click=length(find(Button==3));
                temp_table_clicks.left_click=length(find(Button==1));
                temp_table_clicks.right_click=length(find(Button==3));
               
            catch
                temp_table.left_click=0;
                temp_table.right_click=0;
                temp_table.Flag=0;
                temp_table_clicks.Flag=temp_table.Flag;
                temp_table_clicks.left_click=0;
                temp_table_clicks.right_click=0;
                temp_table_clicks.x{1} = [];
                temp_table_clicks.y{1} = [];
                temp_table_clicks.Button{1} = [];
            end

            
            data_sectors=[data_sectors; temp_table];
            data_clicks=[data_clicks; temp_table_clicks];
                    
            clear temp_table
            clear x y Button


        end
        
        % Save an image with the clicks marked (in cyan) and the radii on
        % top of the GFP picture
        %NOTE: Might want to put these in a 
        saveas(gcf,strcat('Image_Output/',strcat(data.filename{1}(1:end-4),'.png')))
    end
end

%Reorder to get the variables in a nicer order
data_clicks=movevars(data_clicks,'x','Before','y');

save data_sectors_only data_sectors
save data_clicks_all data_clicks
writetable(data_sectors, 'sector_counts.xlsx')
