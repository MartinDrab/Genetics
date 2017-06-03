program VCFGUI;

uses
  Vcl.Forms,
  MainForm in 'MainForm.pas' {MainFrm},
  VCFRecord in 'VCFRecord.pas',
  StatisticsFrm in 'StatisticsFrm.pas' {StatisticsForm};

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TMainFrm, MainFrm);
  Application.Run;
end.
