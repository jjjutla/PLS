% Define channel gains constants
alphasi=1;
alphaid=1;
alphaie=1;

% Initialize an array to store intercept probabilities for Direct Transmission (Pdt), P-AFbORS (Paf), and P-DFbORS (Pdf)
Pdt = zeros(1, 21);
Paf = zeros(1, 21);
Pdf = zeros(1, 21);

% Loop through MER values from -5 to 15
for MER=-5:15
    % Convert MER from dB to linear scale
    lamdade=10^(MER/10);
    
    % Calculate intercept probability for Direct Transmission
    Pdt(MER+6)=1/(1+lamdade);
    
    % Calculate intercept probability for P-AFbORS
    Paf(MER+6)=alphaie/(alphaie+alphaid*lamdade);
    
    % Calculate intercept probability for P-DFbORS
    Pdf(MER+6)=(alphaid+alphasi)/(alphaid+alphasi+alphasi*alphaid*alphaie^-1*lamdade);
end

% Create a plot with semilogarithmic scale on y-axis
figure;

% Define MER values on the x-axis
MER=-5:15

% Plot Direct Transmission intercept probability
x=semilogy(MER,Pdt,'-^');
x.LineWidth=2;
hold on;

% Plot P-AFbORS intercept probability with square exponent
xx=semilogy(MER,Paf.^2,'-x');
xx.LineWidth=2;

% Plot P-DFbORS intercept probability with square exponent
xxx=semilogy(MER,Pdf.^2,'-d');
xxx.LineWidth=2;

% Plot P-AFbORS intercept probability with fourth exponent
xxxx=semilogy(MER,Paf.^4,'-s');
xxxx.LineWidth=2;

% Plot P-DFbORS intercept probability with fourth exponent
xxxxx=semilogy(MER,Pdf.^4,'--');
xxxxx.LineWidth=2;
xxxxx.Color = [0 0 0];

% Turn off the hold for the plot
hold off;
grid on;

% Set plot title, x-axis label, and y-axis label
legend('Direct Transmission', 'P-AFbORS (2 Relays)', 'P-DFbORS (2 relays)', 'P-AFbORS (4 Relays)', 'P-DFbORS (4 Relays)', 'Location', 'best');

ylabel('\fontsize{10}Intercept Probability');
xlabel('\fontsize{10}Main-to-Eavesdropper Ratio MER (dB)');
