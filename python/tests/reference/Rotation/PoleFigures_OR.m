%% Import Script for EBSD Data
%
% Use MTEX
clear ; clear all

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS_bcc = {... 
  crystalSymmetry('m-3m', [2.8665 2.8665 2.8665], 'mineral', 'Iron-alpha', 'color', 'light blue'),...
  crystalSymmetry('m-3m', [1 1 1], 'color', 'light blue')};

CS_fcc = {... 
  crystalSymmetry('m-3m', [3.662 3.662 3.662], 'mineral', 'Iron', 'color', 'light blue'),...
  crystalSymmetry('m-3m', [1 1 1], 'color', 'light blue')};

% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');

%% path to files
pname = 'L:\f.gallardo\DAMASK\python\tests\reference\Rotation'; % has to be changed

% which files to be imported
fname1 = [pname '\bcc_Bain.txt']; fname2 = [pname '\bcc_GT.txt']; fname3 = [pname '\bcc_GT_prime.txt'];
fname4 = [pname '\bcc_KS.txt']; fname5 = [pname '\bcc_NW.txt']; fname6 = [pname '\bcc_Pitsch.txt'];
fname7 = [pname '\fcc_Bain.txt']; fname8 = [pname '\fcc_GT.txt']; fname9 = [pname '\fcc_GT_prime.txt'];
fname10 = [pname '\fcc_KS.txt']; fname11 = [pname '\fcc_NW.txt']; fname12 = [pname '\fcc_Pitsch.txt'];


%% Import the Data

% create an EBSD variable containing the data
ebsd1 = loadEBSD(fname1,CS_bcc,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', 'Passive Rotation');
ebsd2 = loadEBSD(fname2,CS_bcc,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', 'Passive Rotation');
ebsd3 = loadEBSD(fname3,CS_bcc,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', 'Passive Rotation');
ebsd4 = loadEBSD(fname4,CS_bcc,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', 'Passive Rotation');
ebsd5 = loadEBSD(fname5,CS_bcc,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', 'Passive Rotation');
ebsd6 = loadEBSD(fname6,CS_bcc,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', 'Passive Rotation');

ebsd7 = loadEBSD(fname7,CS_fcc,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', 'Passive Rotation');
ebsd8 = loadEBSD(fname8,CS_fcc,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', 'Passive Rotation');
ebsd9 = loadEBSD(fname9,CS_fcc,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', 'Passive Rotation');
ebsd10 = loadEBSD(fname10,CS_fcc,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', 'Passive Rotation');
ebsd11 = loadEBSD(fname11,CS_fcc,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', 'Passive Rotation');
ebsd12 = loadEBSD(fname12,CS_fcc,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', 'Passive Rotation');

%% Plot Data 1stpart_bcc
h1 = [Miller(1,0,0,ebsd1.CS),Miller(1,1,0,ebsd1.CS),Miller(1,1,1,ebsd1.CS)]; % 3 pole figures
plotPDF(ebsd1.orientations,h1,'MarkerSize',5,'MarkerColor','r','MarkerEdgeColor','r','DisplayName','BCC-Bain')
hold on
plotPDF(ebsd2.orientations,h1,'MarkerSize',5,'MarkerColor','w','MarkerEdgeColor','b','DisplayName','BCC-GT')
plotPDF(ebsd3.orientations,h1,'MarkerSize',5,'MarkerColor','w','MarkerEdgeColor','g','DisplayName','BCC-GT_Prime')
legend('show','location','southoutside') 
cd 'L:\f.gallardo\DAMASK\python\tests\reference\Rotation'; % has to be changed 
orient('landscape')
print('-bestfit','1_BCC.pdf','-dpdf')

%% Plot Data 2nd part_bcc
close
plotPDF(ebsd4.orientations,h1,'MarkerSize',5,'MarkerColor','r','MarkerEdgeColor','w','DisplayName','BCC-KS')
hold on
plotPDF(ebsd5.orientations,h1,'MarkerSize',5,'MarkerColor','w','MarkerEdgeColor','m','DisplayName','BCC-NW')
plotPDF(ebsd6.orientations,h1,'MarkerSize',5,'MarkerColor','y','MarkerEdgeColor','w','DisplayName','BCC-Pitsch')
legend('show','location','southoutside')
print('-bestfit','2_BCC.pdf','-dpdf')

%% Plot Data 1stpart_fcc
close
h2 = [Miller(1,0,0,ebsd7.CS),Miller(1,1,0,ebsd7.CS),Miller(1,1,1,ebsd7.CS)]; % 3 pole figures
plotPDF(ebsd7.orientations,h2,'MarkerSize',5,'MarkerColor','r','MarkerEdgeColor','r','DisplayName','FCC-Bain')
hold on
plotPDF(ebsd8.orientations,h2,'MarkerSize',5,'MarkerColor','w','MarkerEdgeColor','b','DisplayName','FCC-GT')
plotPDF(ebsd9.orientations,h2,'MarkerSize',5,'MarkerColor','w','MarkerEdgeColor','g','DisplayName','FCC-GT_Prime')
legend('show','location','southoutside')
print('-bestfit','1_FCC.pdf','-dpdf')

%% Plot Data 2nd part_bcc
close
plotPDF(ebsd10.orientations,h2,'MarkerSize',5,'MarkerColor','r','MarkerEdgeColor','w','DisplayName','FCC-KS')
hold on
plotPDF(ebsd11.orientations,h2,'MarkerSize',5,'MarkerColor','w','MarkerEdgeColor','m','DisplayName','FCC-NW')
plotPDF(ebsd12.orientations,h2,'MarkerSize',5,'MarkerColor','y','MarkerEdgeColor','w','DisplayName','FCC-Pitsch')
legend('show','location','southoutside')
print('-bestfit','2_FCC.pdf','-dpdf')
close

