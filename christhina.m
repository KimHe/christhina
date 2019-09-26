function christhin
% ================================================================== %
% CHRISTHINA: Chromatography RIser THIN Analysis
% ================================================================== %


    j = 0;

    choice = 'No';

    [fileName, pathName, ~] = uigetfile({'*.jpg','JPG (*.jpg)'; ...
        '*.jpeg','JPEG (*.jpeg)'; '*.bmp','BMP Bitmaps'; ...
        '*.tiff','TIFF (*.tiff)'; '*.tif','TIF (*.tif)'; ...
        '*.png','PNG (*.png)'; ...
        '*.*',  'All files (*.*)'}, 'Select a image', 'MultiSelect', 'off');

    [~, folderName, ~] = fileparts(fileName);

    cd(pathName);
    mkdir(folderName);

    figure(1);
    imgRaw = rgb2gray(imread(fileName));
    imshow(imgRaw);

    opts.Interpreter = 'latex';
    opts.Default = 'OK';

    while strcmp(choice,'No') == 1

        rotateAngle = inputdlg({'Input the angle ($$\theta$$) to rotate the image'}, ...
            'Theta Value', [1, 60], {'0'}, opts);

        if str2double(rotateAngle{1}) ~= 0
            imgTLC = imrotate(imgRaw, str2double(cell2mat(rotateAngle)));
        else
            imgTLC = imgRaw;
        end

        imshow(imgTLC);

        choice = questdlg('Is the rotation angle correct?', 'Theta Value', 'OK', 'No', opts);

    end

    cd(folderName);

    imwrite(imgTLC, 'TLCgs.bmp', 'bmp')


    while strcmp(choice, 'OK') == 1

        j = j + 1;

        imgCrop = imcrop(imgTLC);

        imshow(imgCrop);

        [row, ~] = size(imgCrop);
        [~, seedFront] = ginput(2);
        close(gcf);

        if seedFront(1) <= 1; seedFront(1) = 1; end
        if seedFront(2) >= row; seedFront(2) = row-1; end

        seedFront = floor(seedFront);

        imwrite(imgCrop, ['TLCRun' num2str(j) '.bmp'], 'bmp');

        imgCropInv = double( imcomplement(imgCrop) );

        pixels = 1:length(imgCropInv);

        grayScale = mean(imgCropInv.');

        figure(2);
        box on;
        grid on;
        hold all;
        plot(pixels, grayScale, 'k', 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', 'w');
        set(2, 'Color', 'w');
        set(gca, 'FontSize', 20);
        xlabel('Pixels', 'FontSize', 20);
        ylabel('Mean Inverse Grayscale', 'FontSize', 20);
        axis([0 max(pixels) 0 max(grayScale)+10]);

        questdlg('Please select even number of points','Tips');

        [peakGrid, ~] = ginput;

        numPeakGrid = length(peakGrid);
        if mod(numPeakGrid, 2) ~= 0
            error('PlseedFrontYease select even number of points');
        else
            peakNum = numPeakGrid/2;
        end


        areas = zeros(1, peakNum);
        peakIdx = zeros(1, peakNum);
        retentVal = zeros(1, peakNum);

        for i = 1:peakNum

            k = 1 + 2 * (i - 1);

            range = find( pixels < peakGrid(k+1) & pixels > peakGrid(k) );

            x1 = pixels(range);

            q = grayScale(range);

            grad = ((q(1) - q(end)) / (peakGrid(k) - peakGrid(k+1)));

            ord = q(end) - (grad * peakGrid(k+1));

            y1 = grad * x1 + ord;

            areas(i) = trapz(x1, q) - trapz(x1, y1);

            [yMax, idx] = max(grayScale(range));

            xMax = range(idx);

            text(xMax, yMax+5, num2str(i), 'FontSize', 20);

            peakIdx(i) = i;

            retentVal(i) = abs( (max(seedFront) - xMax) / ( max(seedFront) - min(seedFront) ) );

        end


        chroma = [pixels, grayScale];

        str = sprintf(['chromatogram' num2str(j) '.dat']);
        save(str, '-ascii', '-append', '-double', '-tabs');

        print('-f2','-dpdf','chromatogram');

        movefile('chromatogram.pdf', ['chromatogram' num2str(j) '.pdf']);

        close all

        areaTotal = sum(areas);

        areaRatio = (areas / areaTotal) * 100;

        output = [peakIdx; areaRatio; retentVal]';

        %msgbox(num2str(R), ['Report ' num2str(j)],'none')


        TM = round(clock);

        comment = inputdlg('Input a comment','Notes');
        mstxt = ['\n\n' 'Comment: ' cell2mat(comment) '\n'];
        mstxt = [mstxt 'Date: ' num2str(date) '\n'];
        mstxt = [mstxt 'Time: ' num2str(TM(4)) ':' num2str(TM(5)) '\n'];
        mstxt = [mstxt 'Chosen figure: ' fileName '\n' ] ;
        mstxt = [mstxt 'Processed figure: TLCgs.bmp \n'];
        mstxt = [mstxt 'Run: TLCRun' num2str(j) '.bmp \n'];
        mstxt = [mstxt 'Converted to chromatogram' num2str(j) '.pdf \n'];
        mstxt = [mstxt 'Pick\tArea \t\tRf \n' ];

        [r, c] = size(output);
        for k = 1:r
            for l = 1:c
                mstxt = [mstxt num2str(output(k, l)) '\t\t'];
            end
            mstxt = [mstxt '\n'];
        end

        fid = fopen(['output' num2str(j) '.dat'], 'a');
        fprintf(fid,['\n ' mstxt]);
        fclose(fid);

        choice = questdlg('Do you want to realize another run?', ...
            'TLC analysis','OK','No',opts);

    end

    cd('..')

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
