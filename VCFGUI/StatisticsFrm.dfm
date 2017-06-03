object StatisticsForm: TStatisticsForm
  Left = 0
  Top = 0
  Caption = 'Statistics'
  ClientHeight = 348
  ClientWidth = 472
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  Position = poScreenCenter
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object StatisticsPageControl: TPageControl
    Left = 0
    Top = 0
    Width = 472
    Height = 348
    ActivePage = IncidentLengthsTabSheet
    Align = alClient
    TabOrder = 0
    ExplicitWidth = 383
    object IncidentLengthsTabSheet: TTabSheet
      Caption = 'Incident lengths'
      ExplicitLeft = 0
      ExplicitTop = 28
      ExplicitWidth = 375
      object IncidentLengthsListView: TListView
        Left = 0
        Top = 0
        Width = 464
        Height = 320
        Align = alClient
        Columns = <
          item
            Caption = 'Length'
          end
          item
            AutoSize = True
            Caption = 'TPs'
          end
          item
            AutoSize = True
            Caption = 'FNs'
          end
          item
            AutoSize = True
            Caption = 'FPs'
          end>
        ReadOnly = True
        RowSelect = True
        SortType = stBoth
        TabOrder = 0
        ViewStyle = vsReport
        OnCompare = IncidentLengthsListViewCompare
        ExplicitLeft = 24
        ExplicitWidth = 440
      end
    end
  end
end
