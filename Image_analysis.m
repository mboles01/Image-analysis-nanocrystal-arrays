
clear all
close all 

% LOAD IMAGE
nameroot = 'Example_image_2';
eval(['image1 = imread(''' nameroot '.jpg'');']);


% INVERT IMAGE
image2 = imcomplement(image1);
image3 = image2(1:4095,1:4095); 
image4 = imadjust(image3);
image5 = double(image4);


% SCALE IMAGE
% nanometers per pixel, measured using microscopy image scale bar
scale = 20/266;

% BANDPASS FILTER

% Implements a real-space bandpass filter which suppresses pixel noise and 
% long-wavelength image variations while retaining information of 
% characteristic particle size.
% Credit: David Grier (NYU Physics) and Eric Dufresne (Yale Physics)

% out = bpass(image,lnoise,lobject)

% INPUTS:
% image:  The two-dimensional array to be filtered.
% lnoise: Characteristic length scale of noise in pixels. Additive noise 
% averaged over this length should vanish. May assume any positive 
% floating value.
% lobject: A length in pixels somewhat larger than a typical object. 
% Must be an odd valued integer.

lobject = [35]; 
lnoise = [7];

for th = lobject
    for sz = lnoise
        b = bpass(image5,sz,th);
        
%         figure; imagesc(b); set(gcf,'windowstyle','docked','color','w'); colorbar;
%         title(['lnoise =' num2str(qq) ' | lobject = ' num2str(pp) ]);
%         set(gca,'xlim',[1000,2000],'ylim',[1000,2000]); daspect([1,1,1]);
%         
    end
end


% FIND PEAKS
% Credit: Eric Dufresne
% Modified by Michael Boles, March 2014 and February 2018

% out = pkfnd(image,threshold,diameter)
% Provides guess of particle centers by identifying local maxima in
% pixel intensity 
% image: bandpass-filtered (darkfield) image
% threshold: the minimum brightness of a pixel that might be local maximum
% diameter: set to a value slightly larger than particle diameter -- if
% multiple peaks are found within diameter/2, it will keep only the brightest.

% out = cntrd(image,max,diameter)
% Refines particle center guess of pkfnd by taking weighted average of 
% bright spots within particle
% max: locations of local maxima to pixel-level accuracy from pkfnd
% diameter: diamter of window over which to calculate the centroid.
% out: an n-by-4 matrix containing (x,y)-coordinates of particle
% centroids and their intensities

threshold = [30]; 
diameter = [68]; 

for th = threshold
    for diam = diameter
        
        pk = pkfnd(b,th,diam);
        lattice = cntrd(b,pk,diam+1);

        I = lattice(:,3);
        tempidx = find(I > 0);
        lattice = lattice(tempidx,:);

        x = lattice(:,1);
        y = lattice(:,2);
        I = lattice(:,3);
        rsq = lattice(:,4);
        
    end
end


% % Check particle centroid assignment
%         figure; subplot(1,2,1); imagesc(image1); hold on; 
%         hhh = plot(x,y,'or'); set(hhh,'markerfacecolor','r','markersize',3);
%         set(gcf,'windowstyle','docked','color','w');
%         title(['threshold =' num2str(th) ' | diam = ' num2str(diam)  ' | lnoise =' num2str(lnoise) ' | lobject = ' num2str(lobject)  ])
%         daspect([1,1,1]);
% 
%         subplot(1,2,2); hist(I,50); set(get(gca,'child'),'FaceColor',[0.25 0.25 0.25],'EdgeColor','none')
%         title(['threshold =' num2str(th) ' | diam = ' num2str(diam)  ' | lnoise =' num2str(lnoise) ' | lobject = ' num2str(lobject)  ])
%         set(gcf,'windowstyle','docked','color','w'); xlabel('Centroid intensity (a.u.)'); ylabel('Number of particles');

 
% GENERATE INTERPARTICLE BONDS AND MEASURE BOND LENGTHS

% out = bond_analysis(lattice,scale,minlength,maxlength)
% triangulates the set of particle centroids given in lattice, excludes
% spurious bonds by bandpass filtering with minlength & maxlength, which
% are translated from pixels to nanometers with scale

% output is bonds_filtered, a n-by-7 matrix containing:
% (particle no.) (x-coord.) (y-coord.) (bonded to particle no.) (x-coord.) (y-coord.) (bond length in nanometers)

% Written by Michael Boles, University of Chicago, March 2014, updated
% February 2018

minlength = [3];
maxlength = [7];

out = bond_analysis(lattice, minlength, maxlength,scale);
bonds_unfiltered = out{1};
bonds_filtered = out{2};

% % Check bond length distribution for cutoff determination and resulting bond network
% subplot(1,2,1); hist(bonds_unfiltered(:,7),1500); hold on; axis square; set(gca,'xlim',[minlength-1 maxlength+2]); 
% xlabel('Interparticle separation (nm)'); ylabel('Counts'); set(get(gca,'child'),'FaceColor',[0.25 0.25 0.25],'EdgeColor','none');
% subplot(1,2,2); imagesc(image1); daspect([1,1,1]); hold on; 
% for i = 1:length(bonds_filtered)
%     plot([bonds_filtered(i,2),bonds_filtered(i,5)],[bonds_filtered(i,3),bonds_filtered(i,6)],'color','r','linewidth',1);
% end

% COUNT PARTICLE COORDINATION NUMBERS 

% out = coordination_analysis(image5,bonds_filtered,border,lattice);

% Takes bond network contained in bonds_filtered and determines
% coordination state of each particle, excluding particles that are within
% (border) nanometers of the image edge

% Generates a n-by-9 array of containing all the information for particles
% in contact with one another:
% (number) (x-coord.) (y-coord.) (coord. no.) (number) (x-coord.) (y-coord.) (coord. no.) (interparticle distance)

% Also returns n-by-9 arrays of the same data sorted by particle
% coordination endpoints: bonds between 6-fold and 3-fold coordinated particles
% (bonds63), between 6-fold and 4-fold coordinated particles (bonds64), and so on.

% Also returns n-by-3 arrays describing positions of particles with coordination
% numbers 2, 3, 4, 5, and 6 (CN2, CN3, CN4, CN5, and CN6, respectively) in
% the form:  (particle no.) (x-coord) (y-coord)

% Written by Michael Boles, March 2014, updated February 2018

border = [15]; 
out = coordination_analysis(image5,bonds_filtered,border,scale,lattice);
image_analysis_complete = out{10};
CN2 = out{1}; CN3 = out{2}; CN4 = out{3}; CN5 = out{4}; CN6 = out{5};
bonds63 = out{6}; bonds64 = out{7}; bonds65 = out{8}; bonds66 = out{9};

% % plot particles with a given coordination no.
% figure; imagesc(image1); daspect([1,1,1]); hold on; 
% plot(CN2(:,2),CN2(:,3),'.r');
% plot(CN3(:,2),CN3(:,3),'.b');
% plot(CN4(:,2),CN4(:,3),'.g');
% plot(CN5(:,2),CN5(:,3),'.y');
% plot(CN6(:,2),CN6(:,3),'ok');

% plot interparticle bonds separated by coordination number
figure; subplot(1,2,1); imagesc(image1); daspect([1,1,1]); hold on; axis off;
title(nameroot,'Interpreter', 'none');

% bonds between 6- and 3-coord particles
for i = 1:length(bonds63)
    plot([bonds63(i,2), bonds63(i,6)], [bonds63(i,3), bonds63(i,7)],'r');
end
% bonds between 6- and 4-coord particles
for i = 1:length(bonds64)
    plot([bonds64(i,2), bonds64(i,6)], [bonds64(i,3), bonds64(i,7)],'b');
end
% bonds between 6- and 5-coord particles
for i = 1:length(bonds65)
    plot([bonds65(i,2), bonds65(i,6)], [bonds65(i,3), bonds65(i,7)],'g');
end
% bonds between 6- and 6-coord particles
for i = 1:length(bonds66)
    plot([bonds66(i,2), bonds66(i,6)], [bonds66(i,3), bonds66(i,7)],'color',[0.25 0.25 0.25]);
end


% PLOT COORDINATION STATISTICS

set(gcf,'windowstyle','docked','color','w'); 
subplot(4,2,2); hist(bonds66(:,9),50); title('6-6 bonds'); ylabel('Number of bonds'); set(get(gca,'child'),'FaceColor',[0.25 0.25 0.25],'EdgeColor','none'); set(gca,'xlim',[4,7]);
subplot(4,2,4); hist(bonds65(:,9),50); title('6-5 bonds'); ylabel('Number of bonds'); set(get(gca,'child'),'FaceColor','g','EdgeColor','none'); set(gca,'xlim',[4,7]);
subplot(4,2,6); hist(bonds64(:,9),50); title('6-4 bonds'); ylabel('Number of bonds'); set(get(gca,'child'),'FaceColor','b','EdgeColor','none'); set(gca,'xlim',[4,7]);
subplot(4,2,8); hist(bonds63(:,9),50); title('6-3 bonds'); xlabel('Interparticle separation (nm)'); ylabel('Number of bonds'); set(get(gca,'child'),'FaceColor','r','EdgeColor','none'); set(gca,'xlim',[4,7]);

