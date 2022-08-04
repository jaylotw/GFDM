% script on 5/12
% Plot A for slide and thesis
% set default interpreter as latex for the text of figure
% clc; clear; close all;
set(0,'defaultTextInterpreter','latex');

%
% K = 7; 
% M = 8;
N = K*M;
%filterType = 'RC';
%alpha = 0.5;
%[A, g, index_filter] = A_generate(filterType, K, M, alpha);
A_abs = abs(A);
[n, i] = meshgrid(0:N-1);
% F = kron(dftmtx(M)/sqrt(M),eye(K));

% color map
myColormap = [251 251 255;
              236 236 255;
              221 221 255;
              206 206 255;
              185 185 255;
              170 170 255;
              147 147 255;
              125 125 255;
              106 106 255;
              74  74  255;
              40  40  255;
              0   0   227;
              0   0   198;
              0   0   198;
              0   0   147;
              0   0   121];
myColormap = myColormap/255;

% plot A
figure;
surf(n, i, A_abs, 'EdgeColor', [0.6 0.6 0.6]);
colormap(myColormap);
xticks(0:K:N);
yticks(0:K:N);
xlabel('sample index $n$');
ylabel('column index $i$');
xlim([0 N]);ylim([0 N]);
title('$|[\mathbf{A}]_{n,i}|$');
hold on;
% plot g on A
g1 = A(:, 1);
g2 = A(:, 2);
g3 = A(:, 1+K);
plot3(zeros(1,N), 0:N-1, abs(g1), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2);
plot3(ones(1,N), 0:N-1, abs(g2), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
plot3(K*ones(1,N), 0:N-1, abs(g3), 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
set ( gca, 'ydir', 'reverse' );

% plot g
figure;
subplot(3,1,1);
plot(0:N-1, abs(transpose(g1)), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2);
xlim([0 N-1]);
ylabel('$|[\mathbf{A}]_{n,0}|$')
set(gca,'ytick',[])
set(gca,'xtick',[])
box on;

% legend
hold on;
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'--', 'Color', 'black', 'LineWidth', 1);
h(2) = plot(NaN,NaN,':', 'Color', 'black', 'LineWidth', 1);
h(3) = plot(NaN,NaN, 'Color', 'black', 'LineWidth', 1);
hL = legend(h, 'real','imaginary','absolute', 'NumColumns',3);
hL.Location = 'northoutside';

subplot(3,1,2);
hold on;
plot(0:N-1, real(transpose(g2)),'--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
plot(0:N-1, imag(transpose(g2)),':', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
plot(0:N-1, abs(transpose(g2)), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
xlim([0 N-1]);
ylabel('$|[\mathbf{A}]_{n,1}|$')
set(gca,'xtick',[])
set(gca,'ytick',[])
box on;

subplot(3,1,3);
plot(0:N-1, abs(transpose(g3)), 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
xlim([0 N-1]);
ylabel('$|[\mathbf{A}]_{n,K}|$')
set(gca,'ytick',[])
xlabel('sample index $n$');
box on;