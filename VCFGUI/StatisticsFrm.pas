Unit StatisticsFrm;

Interface

Uses
  Winapi.Windows, Winapi.Messages, System.SysUtils,
  System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.ComCtrls,
  VCFRecord, Generics.Collections;

Type
  TIncidentLengthStat = Record
    TPs : Cardinal;
    FNs : Cardinal;
    FPs : Cardinal;
    end;
  PTIncidentLengthStat = ^TIncidentLengthStat;

  TStatisticsForm = Class(TForm)
    StatisticsPageControl: TPageControl;
    IncidentLengthsTabSheet: TTabSheet;
    IncidentLengthsListView: TListView;
    procedure FormCreate(Sender: TObject);
    procedure IncidentLengthsListViewCompare(Sender: TObject; Item1,
      Item2: TListItem; Data: Integer; var Compare: Integer);
  Private
    FIncidents : TList<TVCFIncident>;
    FTruthSet : TList<TVCFRecord>;
    FTestSet : TList<TVCFRecord>;

    FILs : TDictionary<Cardinal, PTIncidentLengthStat>;
    Procedure IncidentLengthStatistics;
  Public
    Constructor Create(AOwner:TComponent; ATruthSet:TList<TVCFRecord>; ATestSet:TList<TVCFRecord>; AIncidents:TList<TVCFIncident>); Reintroduce;
  end;


Implementation

{$R *.DFM}

Constructor TStatisticsForm.Create(AOwner:TComponent; ATruthSet:TList<TVCFRecord>; ATestSet:TList<TVCFRecord>; AIncidents:TList<TVCFIncident>);
begin
FIncidents := AIncidents;
FTruthSet := ATruthSet;
FTestSet := ATestSet;
Inherited Create(AOwner);
end;

Procedure TStatisticsForm.FormCreate(Sender: TObject);
begin
FILs := TDictionary<Cardinal, PTIncidentLengthStat>.Create;
IncidentLengthStatistics;
end;

Procedure TStatisticsForm.IncidentLengthsListViewCompare(Sender: TObject; Item1,
  Item2: TListItem; Data: Integer; var Compare: Integer);
begin
Compare := Integer(Item1.Data) - Integer(Item2.Data);
end;

Procedure TStatisticsForm.IncidentLengthStatistics;
Var
  len : Cardinal;
  i : TVCFIncident;
  p : TPair<Cardinal, PTIncidentLengthStat>;
  il : PTIncidentLengthStat;
  totalTPs : Cardinal;
  totalFNs : Cardinal;
  totalFPs : Cardinal;
  InILs : Boolean;
begin
totalTPs := 0;
totalFPs := 0;
totalFNs := 0;
For i In FIncidents Do
  begin
  If (I.VCFRecord.RecordType = vcfrtInsertion) Or
    (I.VCFRecord.RecordType = vcfrtDeletion) Then
    begin
    Case I.VCFRecord.RecordType Of
      vcfrtInsertion : len := Length(I.VCFRecord.Alt);
      vcfrtDeletion : len := Length(I.VCFRecord.Ref);
      end;

    inILs := FILs.TryGetValue(len, il);
    If Not inILs Then
      begin
      New(il);
      il.TPs := 0;
      il.FNs := 0;
      il.FPs := 0;
      end;

    Case I.IncidentType Of
      vcfitTruePositive : begin
        Inc(il.TPs);
        Inc(totalTPs);
        end;
      vcfitFalseNegative : begin
        Inc(il.FNs);
        Inc(totalFNs);
        end;
      vcfitFalsePositive : begin
        Inc(il.FPs);
        Inc(totalFPs);
        end;
      end;

    If Not inILs Then
      FILs.Add(len, il);
    end;
  end;

With IncidentLengthsListVIew Do
  begin
  Clear;
  Items.BeginUpdate;
  For p In FILs Do
    begin
    With Items.Add Do
      begin
      Caption := Format('%u', [p.Key]);
      SubItems.Add(Format('%u (%f %%)', [p.Value.TPs, p.Value.TPs / totalTPs * 100]));
      SubItems.Add(Format('%u (%f %%)', [p.Value.FNs, p.Value.FNs / totalFNs * 100]));
      SubItems.Add(Format('%u (%f %%)', [p.Value.FPs, p.Value.FPs / totalFPs * 100]));
      Data := Pointer(p.Key);
      end;
    end;

  Items.EndUpdate;
  end;
end;


End.
