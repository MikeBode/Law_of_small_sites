function Rare_species_figure()

% ================================================================================================
% This function analyses the data and produces the outputs described in:
%   The law of small sites and its implications for conservation
%   Mahood et al. (submitted).
% 
% This dataset is sourced from Dryad at:
% Bovendorp RS, Brum FT, McCleery RA, Baiser B, Loyola R, Cianciaruso MV, Galetti M (2018) 
% Data from: Defaunation and fragmentation erode small mammal diversity dimensions in tropical forests. 
% Dryad Digital Repository. https://doi.org/10.5061/dryad.p7h807v

% A description of the dataset can be found in:
% Bovendorp RS, Brum FT, McCleery RA, Baiser B, Loyola R, Cianciaruso MV, Galetti M (2018) 
% Defaunation and fragmentation erode small mammal diversity dimensions in tropical forests. 
% Ecography, online in advance of print.
% ================================================================================================

% Import the empirical data from an Excel spreadsheet. 
[D,T] = xlsread('Empiricaldata.xlsx');

% The area of each patch is in the final column
Area = D(:,end);

% The species-site matrix is the rest of the datasheet
SSM = D(:,1:end-1);

% Extract the number of patches (rows) and species (columns)
NumPatches = size(SSM,1)
NumSpecies = size(SSM,2)

% Randomly re-sample the dataset
for s = 1:NumSpecies
   
   % For each species, re-allocate the abundance in each patch to 
   % a different patch at random.
   SSM_bs(:,s) = randsample(SSM(:,s),NumPatches);
end

% Create the outputs for both the empirical data and the randomly
% resampled version.
figure(1), clf
sub_MakePanels(SSM,Area,1)
sub_MakePanels(SSM_bs,Area,2)

function sub_MakePanels(SSM,Area,TE)
% This subroutine creates the output figures. Inputs are 
%     SSM - the species site matrix
%     Area - a column vector of site areas
%     TE - this variable indicates whether we're graphing the empirical data (1) or the simulated data (2)

% What size font for the figure text and labels
FS = 16;

% Re-extract the number of species and patches
NumSpecies = size(SSM,2)
NumPatches = size(SSM,1)

% Define a vector that indicates how many patches each species is found on
Dist = sum(SSM>0);

% For each species, calculate the average size of the patches it occupies.
for s = 1:NumSpecies
   F = find(SSM(:,s)>0);
   AvPatchSize(s) = mean(Area(F));
end

if TE == 1
   % Create a histogram of the different fragment sizes
   subplot(4,2,1); box on
   X = linspace(min(Area),max(Area),50);
   H = histc(Area,X); H = H./sum(H);
   bb = bar(X,H,0.8); set(bb,'facecolor',hex2rgb('#089404'));
   xlabel('Fragment size','fontsize',FS)
   ylabel('Frequency','fontsize',FS); axis tight
   xlim([-1.2e4 2.1e5]); ylim([0 1]), box on
   
   % Create a histogram of the different species' abundances
   subplot(4,2,2); box on
   X = linspace(min(Dist),max(Dist),20);
   H = histc(Dist,X); H = H./sum(H);
   bb = bar(X,H,0.7); set(bb,'facecolor',hex2rgb('#089404'));
   xlabel('Species range (# fragments)','fontsize',FS)
   ylabel('Frequency','fontsize',FS); axis tight
   xlim([-5 155]); ylim([0 0.6]), box on
end

% Plot rare species fragment sizes versus all fragment sizes
if TE == 1; subplot(4,2,3); hold on % If we're plotting the empirical data
else; subplot(4,2,4); hold on % If we're plotting the randomised data
end

% Define a "rare" species. 
% For example, if RARE = 1, then a rare species is only found in one patch
RARE = 1; 

% Log transform the habitat area
LA = log(Area);

% Indentify the species that are "rare"
F = find(sum(SSM>0)<=RARE);

% What's the area of the patches that these species inhabit?
LA_endemics = [];
for f = 1:length(F)
   Ff = find(SSM(:,F(f))>0);
   LA_endemics = [LA_endemics; LA(Ff)];
end

% Histogram of the area of all patches, and the area of patches with rare species
X = linspace(min(LA),max(LA),20);
Ha = histc(LA,X);
Hr = histc(LA_endemics,X);
B = bar(X,Ha./sum(Ha),1); set(B,'facecolor','none','linewidth',1)
B = bar(X,Hr./sum(Hr),0.5); set(B,'facecolor','r','edgecolor','none')
xlim([-3 14]); box off
l = legend('All fragments','Rare species');
set(l,'box','off','location','northwest')
ylim([0 0.4])

% Perform a K-S test to see if the patches with rare species are a random sample from all patches
[h,p] = kstest2(LA,LA_endemics);
disp(['p-value = ' num2str(p,3)])

% Plot the relative density vs fragment size graph
if TE == 1; subplot(4,2,5); hold on; box on % If we're plotting the empirical data
else;       subplot(4,2,6); hold on; box on % If we're plotting the randomised data
end


% Convert species' abundances into an integer color-scheme from 1-grads
grads = 25;
CL = flipud(jet(grads));
RawAbundance = sum(SSM);
Abundance = log(sum(SSM));
Abundance = sum(SSM);
Abundance = Abundance - min(Abundance);
Abundance = Abundance./max(Abundance);
Abundance = 1 + ceil((grads-1)*Abundance);

% Calculate global abundance and global density for each species
GlobalAbundance = sum(SSM);
GlobalDensity = GlobalAbundance./sum(Area);

% Go through all sites and all species, and create the data for a scatter plot
xx = []; yy = []; CLs = [];
for st = 1:size(SSM,1)
   for spp = 1:size(SSM,2)
      
      % This is the density of the species in this site, compared to that species' global density
      yy = [yy SSM(st,spp)./Area(st)./GlobalDensity(spp)];
      
      % This the log of the area of that same site
      xx = [xx log(Area(st))];
      
      % This is the color of the marker (based on the species global abundance
      CLs = [CLs; CL(Abundance(spp),:)];
   end
end

% Scatter plot the information above
scatter(xx,log(yy),30,CLs,'.'); axis tight
xlabel('Patch area (ln)','fontsize',FS)
ylabel('Local density (ln)','fontsize',FS)
XL = xlim; YL = ylim; xlim(XL+[-1 1]); ylim(YL+[-1 1])

% Plot the cumulative rare species found in the smallest % of sites

% Define and find the "rare" species
RareDef = 1;
F_rare = find(sum(SSM>0)<=RareDef);
NumRareSp = length(F_rare);

% Go through the sites by percentage of the mean site
T = 1000; C = round(T*0.08);
for i = 1:C
   
   % How many sites are smaller than this percentile?
   F_small = find(Area < i*mean(Area)/T);
   
   % How many unique species are found in those small sites?
   NumRare(i) = sum(sum(SSM(F_small,F_rare))>0);
end

if TE == 1; subplot(4,2,7); hold on; box on % If we're plotting the empirical data
else;       subplot(4,2,8); hold on; box on % If we're plotting the randomised data
end

% Create a color for the plot (red-ish)
CLs = get(gca,'colororder');
pp = patch(100*[0 1:C C 0]./T,100*[0 NumRare./NumRareSp 0 0],CLs(7,:));
set(pp,'facealpha',0.3,'linewidth',2,'edgecolor',CLs(7,:))
set(gca,'fontsize',FS-4,'ytick',[0:25:100])
xlabel('Patch area (% mean)','fontsize',FS)
ylabel('% rare species','fontsize',FS)
ylim([0 100])

