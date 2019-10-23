function christhina
% ================================================================== %
% CHRISTHINA: CHromatography RISer THIN Analysis
% ================================================================== %


    j = 0;

    iterative = 'Yes';

    imgTLC = imgPreProcess;

    while strcmp(iterative, 'Yes') == 1

        j = j + 1;

        imgCrop = imcrop(imgTLC);

        imshow(imgCrop);

        [~, seedFront] = ginput(2);

        close(gcf);

        if seedFront(1) <= 1; seedFront(1) = 1; end

        if seedFront(2) >= size(imgCrop)
            seedFront(2) = size(imgCrop) - 1;
        end

        seedFront = floor(seedFront);

        imwrite(imgCrop, ['TLCRun' num2str(j) '.bmp'], 'bmp');

        imgCropInv = double( imcomplement(imgCrop) );

        pixels = 1:length(imgCropInv);

        grayScale = mean(imgCropInv.');

        [peakNum, peakGrid] = imgPlotPick(pixels, grayScale);


        areas = zeros(1, peakNum);
        peakIdx = zeros(1, peakNum);
        retentVal = zeros(1, peakNum);

        for i = 1:peakNum

            k = 1 + 2 * (i - 1);

            range = find( pixels < peakGrid(k+1) & pixels > peakGrid(k) );

            x = pixels(range);

            y = grayScale(range);

            grad = ( y(1) - y(end) ) / ( peakGrid(k) - peakGrid(k+1) );

            ord = y(end) - grad * peakGrid(k+1);

            y1 = grad * x + ord;

            %areas(i) = trapz(x, y) - trapz(x, y1);
            areas(i) = trapz(x, y);

            [yMax, idx] = max(grayScale(range));

            xMax = range(idx);

            text(xMax, yMax+5, num2str(i), 'FontSize', 20);

            peakIdx(i) = i;

            retentVal(i) = abs( (max(seedFront) - xMax) / ( max(seedFront) - min(seedFront) ) );

        end

        areaTotal = sum(areas);

        areaRatio = (areas / areaTotal) * 100;

        output = [peakIdx; areaRatio; retentVal]';

        dataWrite(j, pixels, grayScale, output);

        iterative = questdlg('Do you want to run repeatedly?', ...
            'TLC analysis', 'Yes', 'No', 'Yes');

    end

    cd('..')

end


function imgTLC = imgPreProcess
% ================================================================== %
% Read rotate and show the loaded image
% ================================================================== %


    global fileName;

    rotate = 'No';

    fileName = uigetfile({'*.jpg','JPG (*.jpg)'; '*.jpeg','JPEG (*.jpeg)';
    '*.bmp','BMP Bitmaps'; '*.tiff','TIFF (*.tiff)'; '*.tif','TIF (*.tif)'; ...
        '*.png','PNG (*.png)'; '*.*', 'All files (*.*)'}, ...
        'Select a image', 'MultiSelect', 'off');

    [~, folderName, ~] = fileparts(fileName);

    if exist(folderName, 'dir') ~= 7
        mkdir(folderName);
    end

    figure(1);

    imgRaw = rgb2gray( imread(fileName) );

    imshow(imgRaw);

    opts.Interpreter = 'latex';
    opts.Default = 'Yes';

    while strcmp(rotate, 'No') == 1

        rotateAngle = inputdlg({'Input the angle ($$\theta$$) to rotate the image'}, ...
            'Theta Value', [1, 60], {'0'}, opts);

        if str2double(rotateAngle) ~= 0
            imgTLC = imrotate(imgRaw, str2double(cell2mat(rotateAngle)));
        else
            imgTLC = imgRaw;
        end

        imshow(imgTLC);

        rotate = questdlg('Is the rotation angle correct?', 'Theta Value', 'Yes', 'No', opts);

    end

    imgTLC = imgRaw;

    cd(folderName);

    imwrite(imgTLC, 'TLCgs.bmp', 'bmp')

end


function [peakNum, peakGrid] = imgPlotPick(pixels, grayScale)
% ================================================================== %
% Plot chromatogram and pick up the peaks to be analyzed
% ================================================================== %


    figure(2);

    box on;
    grid on;
    hold all;

    plot(pixels, grayScale, 'k', 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', 'w');

    axis([0 max(pixels) 0 max(grayScale)+10]);

    set(2, 'Color', 'w');
    set(gca, 'FontSize', 20);

    xlabel('Pixels', 'FontSize', 20);
    ylabel('Mean Inverse Grayscale', 'FontSize', 20);


    questdlg('Please select even number of points', 'Tips');

    [peakGrid, ~] = ginput;

    numPeakGrid = length(peakGrid);

    if mod(numPeakGrid, 2) ~= 0
        error('Plsease select even number of points');
    else
        peakNum = numPeakGrid/2;
    end

end


function dataWrite(j, pixels, grayScale, output)
% ================================================================== %
% Write out the chromatogram and log file
% ================================================================== %


    global fileName;

    chroma = [pixels', grayScale'];

    str = sprintf(['chromatogram' num2str(j) '.dat']);
    save(str, 'chroma', '-ascii', '-tabs');

    print('-f2','-dpdf','chromatogram');
    movefile('chromatogram.pdf', ['chromatogram' num2str(j) '.pdf']);

    close all


    TM = round(clock);

    comment = inputdlg('Input a comment','Notes');

    strOut = sprintf('Comment: %s \n Date: %s \n Time: %s : %s \n Chosen figure: %s \n Processed figure TLCgs.bmp \n Run: TLCRun%s.bmp \n Converted to chromatogram%s.pdf \n Pick\tArea \t\tRetention value \n',...
        cell2mat(comment), num2str(date), num2str(TM(4)), num2str(TM(5)), fileName, num2str(j), num2str(j));

    fid = fopen(['output' num2str(j) '.dat'], 'a');
    fprintf(fid,['\n ' strOut]);

    fclose(fid);

    save(sprintf('output%s.dat', num2str(j)), 'output', '-ascii', '-append', '-tabs');

end
% ================================================================== %
%                   Copyright (C) 2013-2019
%                    Maximiliano Barchiesi
%                       Mercedes Sangroni
%                        Carlos Renaudo
%               Pablo Rossi <simap@ing.unrc.edu.ar>
%              Qiao-Le He <qiaole.he@rwth-aachen.de>
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% ================================================================== %
